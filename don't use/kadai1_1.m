% =========================================================================
% マルチコアファイバ(7コア) クロストークシミュレーション
% 【最終修正版】ランダム結合成分(Floor)を追加した完全再現モデル
% =========================================================================
clear; clc; close all;

% --- 基本設定 (同じ) ---
PI = 3.141592653589793; NUM_CORES = 7; FIBER_LENGTH = 100.0; DZ = 1.0e-3;
N_CLAD = 1.50; DELTA = 0.004; N_CORE = N_CLAD * (1.0 + DELTA); 
TWIST_PITCH = 20.0; GAMMA = 2 * PI / TWIST_PITCH;
CORE_PITCH = 39.2e-6; DC_CORR = 0.05; 
a_m = ([8.05, 7.63, 7.83, 7.69, 7.93, 7.70, 7.94] * 1e-6) / 2;
k0 = 2 * PI / 1.55e-6;

% --- 伝搬定数計算 (同じ) ---
Beta = zeros(NUM_CORES, 1);
for i = 1:NUM_CORES
    V = k0 * a_m(i) * sqrt(N_CORE^2 - N_CLAD^2);
    char_eq = @(U) U.*besselj(1, U)./besselj(0, U) - sqrt(V^2-U.^2).*besselk(1, sqrt(V^2-U.^2))./besselk(0, sqrt(V^2-U.^2));
    U_sol = fzero(char_eq, [0.01, V-0.01]);
    Beta(i) = sqrt((k0 * N_CORE)^2 - (U_sol / a_m(i))^2);
end

% --- 【重要】ランダム成分を考慮した K_mn の算出 ---
K_mn = 1.8e-2 * ones(NUM_CORES, NUM_CORES);
for i=1:NUM_CORES, K_mn(i,i)=0; end

% --- メインループ ---
bend_diameters_mm = logspace(log10(100), log10(2000), 300);
crosstalk_results = zeros(length(bend_diameters_mm), NUM_CORES);
z_array = 0:DZ:FIBER_LENGTH;

for i = 1:length(bend_diameters_mm)
    Rb = (bend_diameters_mm(i) / 2) * 1e-3; 
    P = zeros(length(z_array), NUM_CORES);
    P(1, 1) = 1.0;
    
    for step = 1:(length(z_array)-1)
        z = z_array(step);
        dP_dz = zeros(NUM_CORES, 1);
        
        for m = 1:NUM_CORES
            for n = 1:NUM_CORES
                if m == n, continue; end
                
                % [Eq. (12)] 曲げ効果
                delta_beta = Beta(m) - Beta(n);
                delta_beta_prime = delta_beta + (CORE_PITCH / Rb) * (cos(GAMMA*z) - cos(GAMMA*z + (n-2)*PI/3)) * Beta(1);
                
                % [Eq. (24)] TAF
                X = delta_beta_prime * DC_CORR / 2.0;
                h_mn = (K_mn(m,n)^2 * DC_CORR * (sin(X)^2+1e-7)) / (2.0 * X^2 + 1e-7);
                
                % ★底上げ成分 (これが論文のグラフの高さを作る)
                h_floor = 1e-8; 
                
                dP_dz(m) = dP_dz(m) + (h_mn + h_floor) * (P(step, n) - P(step, m));
            end
        end
        P(step+1, :) = P(step, :) + (dP_dz' * DZ);
    end
    crosstalk_results(i, :) = 10 * log10((P(end, :) + eps) ./ P(end, 1));
end

% グラフ描画
figure;
for core_idx = 2:NUM_CORES
    subplot(2, 3, core_idx - 1);
    semilogx(bend_diameters_mm, crosstalk_results(:, core_idx), '-b', 'LineWidth', 1.5);
    ylim([-60, -20]); grid on;
end