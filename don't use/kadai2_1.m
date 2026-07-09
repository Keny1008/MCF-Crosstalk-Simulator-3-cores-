% =========================================================================
% マルチコアファイバ(7コア) クロストークシミュレーション
% 結合電力理論(CPT) + 指数自己相関関数(EAF) [Fig. 5 再現版]
% =========================================================================
clear; clc; close all;

% --- 1. 基本パラメータ設定 ---
PI = 3.141592653589793;          
NUM_CORES = 7;                   
FIBER_LENGTH = 100.0;            
DZ = 1.0e-3;                     
N_CLAD = 1.50;                   
DELTA = 0.004; 
N_CORE = N_CLAD * (1.0 + DELTA); 
TWIST_PITCH = 20.0;              
E0 = 8.85418781e-12;             
M0 = 1.25663706e-6;              
C_SPEED = 1.0 / sqrt(E0 * M0);   
F = 193e12;                      
WAVELENGTH = C_SPEED / F;        
CORE_PITCH = 39.2e-6;            

% 各コアの半径 [m]
CORE_DIAMETERS = [8.05, 7.63, 7.83, 7.69, 7.93, 7.70, 7.94] * 1e-6;
a_m = CORE_DIAMETERS / 2;        

gamma = 2 * PI / TWIST_PITCH;    
k0 = 2 * PI / WAVELENGTH;        

% --- 2. 伝搬定数 Beta と モード結合係数 K_mn の理論計算 ---
Beta = zeros(NUM_CORES, 1);
U_vec = zeros(NUM_CORES, 1);
W_vec = zeros(NUM_CORES, 1);
V_vec = zeros(NUM_CORES, 1);

fprintf('伝搬定数と結合係数を計算中...\n');
for i = 1:NUM_CORES
    V = k0 * a_m(i) * sqrt(N_CORE^2 - N_CLAD^2);
    V_vec(i) = V;
    
    char_eq = @(U) U .* besselj(1, U) ./ besselj(0, U) - ...
                   sqrt(V^2 - U.^2) .* besselk(1, sqrt(V^2 - U.^2)) ./ besselk(0, sqrt(V^2 - U.^2));
               
    options = optimset('Display', 'off');
    U_sol = fzero(char_eq, [0.01, V - 0.01], options);
    W_vec(i) = sqrt(V^2 - U_sol^2);
    U_vec(i) = U_sol;
    
    Beta(i) = sqrt((k0 * N_CORE)^2 - (U_sol / a_m(i))^2);
end

dist_matrix = zeros(NUM_CORES, NUM_CORES);
for m = 1:NUM_CORES
    for n = 1:NUM_CORES
        if m == n, continue; end
        if m == 1 || n == 1
            dist_matrix(m,n) = CORE_PITCH; 
        else
            angle_m = (m - 2) * (PI / 3);
            angle_n = (n - 2) * (PI / 3);
            dist_matrix(m,n) = sqrt(2*CORE_PITCH^2 - 2*CORE_PITCH^2*cos(angle_m - angle_n));
        end
    end
end

K_mn = zeros(NUM_CORES, NUM_CORES);
for m = 1:NUM_CORES
    for n = 1:NUM_CORES
        if m == n, continue; end
        
        a_avg = (a_m(m) + a_m(n)) / 2;
        V_avg = (V_vec(m) + V_vec(n)) / 2;
        U_avg = (U_vec(m) + U_vec(n)) / 2;
        W_avg = (W_vec(m) + W_vec(n)) / 2;
        d = dist_matrix(m,n);
        
        kappa = (sqrt(2 * DELTA) / a_avg) * (U_avg^2 / V_avg^3) * ...
                (besselk(0, W_avg * d / a_avg) / besselk(1, W_avg)^2);
        K_mn(m,n) = kappa;
    end
end

% 理論値と現実のズレの補正係数 (前回と同じく維持)
K_mn = K_mn * 2.2; 

% --- 3. メインループの準備 ---
% 図5を再現するため、4つの相関長を用意
dc_array = [0.01, 0.05, 0.1, 0.5];
line_styles = {':m', '-b', '--g', '-.k'}; % マゼンタ点線, 青実線, 緑破線, 黒一点鎖線

bend_diameters_mm = logspace(log10(100), log10(2000), 300);
crosstalk_results = zeros(length(dc_array), length(bend_diameters_mm), NUM_CORES);
z_array = 0:DZ:FIBER_LENGTH;
num_steps = length(z_array);

fprintf('シミュレーションを開始します... (4パターンの相関長を計算します)\n');

% 相関長(dc)ごとの外側ループ
for dc_idx = 1:length(dc_array)
    DC_CORR = dc_array(dc_idx);
    fprintf('  >> dc = %.2f m の計算中...\n', DC_CORR);
    
    for i = 1:length(bend_diameters_mm)
        Rb = (bend_diameters_mm(i) / 2) * 1e-3; 
        
        P = zeros(num_steps, NUM_CORES);
        P(1, 1) = 1.0;
        
        for step = 1:(num_steps-1)
            z = z_array(step);
            dP_dz = zeros(NUM_CORES, 1);
            
            for m = 1:NUM_CORES
                for n = 1:NUM_CORES
                    if m == n, continue; end
                    
                    % [Eq. (4)] 
                    theta_m = gamma * z + (m - 2) * (PI / 3);
                    theta_n = gamma * z + (n - 2) * (PI / 3);
                    cos_m = 0; if m > 1, cos_m = cos(theta_m); end
                    cos_n = 0; if n > 1, cos_n = cos(theta_n); end
                    
                    % [Eq. (12)] 
                    delta_beta = Beta(m) - Beta(n); 
                    delta_beta_prime = delta_beta + (CORE_PITCH / Rb) * (Beta(m) * cos_m - Beta(n) * cos_n);
                    
                    % --- ★変更箇所: [Eq. (22)] EAFによる電力結合係数 h_mn ---
                    h_mn = (K_mn(m,n)^2 * DC_CORR) / (1.0 + (delta_beta_prime * DC_CORR)^2);
                    
                    % [Eq. (11)] 
                    dP_dz(m) = dP_dz(m) + h_mn * (P(step, n) - P(step, m));
                end
            end
            P(step+1, :) = P(step, :) + (dP_dz' * DZ);
        end
        crosstalk_results(dc_idx, i, :) = 10 * log10((P(end, :) + eps) ./ P(end, 1));
    end
end

fprintf('シミュレーション完了！グラフを描画します。\n');

% --- 4. グラフ描画 ---
figure('Position', [100, 100, 900, 600]); 
for core_idx = 2:NUM_CORES
    subplot(2, 3, core_idx - 1); 
    
    % 相関長ごとに異なる線種で重ね描き
    for dc_idx = 1:length(dc_array)
        semilogx(bend_diameters_mm, squeeze(crosstalk_results(dc_idx, :, core_idx)), ...
                 line_styles{dc_idx}, 'LineWidth', 1.2);
        hold on;
    end
    
    grid on;
    title(sprintf('Core %d', core_idx));
    xlabel('Bending diameter [mm]');
    ylabel('Crosstalk [dB]');
    xlim([100, 2000]);
    ylim([-60, -20]);
    xticks([100, 200, 500, 1000, 2000]); 
end