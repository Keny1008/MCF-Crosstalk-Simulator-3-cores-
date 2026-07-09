% =========================================================================
% マルチコアファイバ(7コア) クロストークの曲げ直径依存性シミュレーション
% 結合電力理論 (CPT) + 三角自己相関関数 (TAF) [論文パラメータ・フィッティング版]
% =========================================================================
clear; clc; close all;

% --- 1. 基本パラメータ設定 ---
PI = 3.141592653589793;          
NUM_CORES = 7;                   
FIBER_LENGTH = 100.0;            % ファイバ長 [m]
DZ = 1.0e-3;                     % 空間ステップ幅 [m]
N_CLAD = 1.50;                   % クラッド屈折率
N_CORE = N_CLAD * (1.0 + 0.003); % コア屈折率 (0.3%)
TWIST_PITCH = 20.0;              % ツイストピッチ [m]
E0 = 8.85418781e-12;             
M0 = 1.25663706e-6;              
C_SPEED = 1.0 / sqrt(E0 * M0);   
F = 193e12;                      % 周波数 [Hz]
WAVELENGTH = C_SPEED / F;        % 波長 [m]
DC_CORR = 0.05;                  % 相関長 [m] 
CORE_PITCH = 39.2e-6;            % コア間距離 [m]

% ツイストの角速度 [rad/m]
gamma = 2 * PI / TWIST_PITCH;

% --- 2. 事前準備（伝搬定数 Beta と モード結合係数 K_mn） ---
k0 = 2 * PI / WAVELENGTH;

% 中心コア(Core 1)の伝搬定数を基準とする
Beta_1 = N_CORE * k0; 

% 論文のFig.9のピーク位置から逆算した各コアと中心コアの伝搬定数差 [rad/m]
% (Core1との差: Core1=0, Core2=800, Core3=410, Core4=690, Core5=270, Core6=670, Core7=250)
DeltaBeta = [0, 800, 410, 690, 270, 670, 250]; 

% 各コアの伝搬定数 (中心コアが一番太いので、外周コアはBetaが小さくなる)
Beta = Beta_1 - DeltaBeta;

% モード結合係数 K_mn の設定
% (論文のベースライン -50dB 付近に合わせるための調整値)
K_mn = 1.8e-2 * ones(NUM_CORES, NUM_CORES); 
for i = 1:NUM_CORES
    K_mn(i,i) = 0; % 自分自身への結合はゼロ
end

% --- 3. 曲げ直径の配列作成と保存用変数の準備 ---
% 100mm から 2000mm まで対数間隔で50点のデータポイントを作成
bend_diameters_mm = logspace(log10(100), log10(2000), 50);

% 計算結果を保存する配列
crosstalk_results = zeros(length(bend_diameters_mm), NUM_CORES);

% Z軸のグリッド
z_array = 0:DZ:FIBER_LENGTH;
num_steps = length(z_array);

% --- 4. メインループ (曲げ直径ごとに計算) ---
fprintf('シミュレーションを開始します... (全50ステップ)\n');

for i = 1:length(bend_diameters_mm)
    % 直径[mm]から半径[m]に変換
    Rb = (bend_diameters_mm(i) / 2) * 1e-3; 
    
    % 電力の初期化 (z=0 で中心コアに電力1を投入)
    P = zeros(num_steps, NUM_CORES);
    P(1, 1) = 1.0;
    
    % Z方向への伝搬シミュレーション
    for step = 1:(num_steps-1)
        z = z_array(step);
        dP_dz = zeros(NUM_CORES, 1);
        
        for m = 1:NUM_CORES
            for n = 1:NUM_CORES
                if m == n, continue; end
                
                % --- [Eq. (4)] ツイストによる各コアの回転角 ---
                theta_m = gamma * z + (m - 2) * (PI / 3);
                theta_n = gamma * z + (n - 2) * (PI / 3);
                
                cos_m = 0; if m > 1, cos_m = cos(theta_m); end
                cos_n = 0; if n > 1, cos_n = cos(theta_n); end
                
                % --- [Eq. (12)] 局所的な伝搬定数差 ---
                delta_beta = Beta(m) - Beta(n); 
                delta_beta_prime = delta_beta + (CORE_PITCH / Rb) * (Beta(m) * cos_m - Beta(n) * cos_n);
                
                % --- [Eq. (24)] TAFによる電力結合係数 h_mn ---
                X = delta_beta_prime * DC_CORR / 2.0;
                if abs(X) < 1e-12
                    h_mn = (K_mn(m,n)^2 * DC_CORR) / 2.0;
                else
                    h_mn = (K_mn(m,n)^2 * DC_CORR * sin(X)^2) / (2.0 * X^2);
                end
                
                % --- [Eq. (11)] 結合電力方程式 ---
                dP_dz(m) = dP_dz(m) + h_mn * (P(step, n) - P(step, m));
            end
        end
        % パワー更新
        P(step+1, :) = P(step, :) + (dP_dz' * DZ);
    end
    
    % 100m地点のクロストーク[dB]を計算して保存
    crosstalk_results(i, :) = 10 * log10((P(end, :) + eps) ./ P(end, 1));
end

fprintf('シミュレーション完了！グラフを描画します。\n');

% --- 5. 論文と同じ構成でのグラフ描画 ---
figure('Position', [100, 100, 900, 600]); 

for core_idx = 2:NUM_CORES
    subplot(2, 3, core_idx - 1); 
    
    % x軸を対数スケールにしてプロット (青色の実線)
    semilogx(bend_diameters_mm, crosstalk_results(:, core_idx), '-b', 'LineWidth', 1.5);
    hold on;
    
    grid on;
    title(sprintf('Core %d', core_idx));
    xlabel('Bending diameter [mm]');
    ylabel('Crosstalk [dB]');
    
    % 表示範囲と目盛りを論文に合わせる
    xlim([100, 2000]);
    ylim([-60, -20]);
    xticks([100, 200, 500, 1000, 2000]); 
end