% =========================================================================
% 4コアMCF クロストーク解析：ファイバ長 100m 実装
% =========================================================================
clear; clc; close all;

% --- 1. パラメータ設定 ---
NUM_CORES = 4;
FIBER_LENGTH = 100.0; % [m] 指定通り
DZ = 0.5;             % 計算ステップ [m]
CORE_PITCH = 39.2e-6; % [m]
DC_CORR = 0.05;       % 相関長 dc = 5cm [cite: 176, 242]
WAVELENGTH = 1.55e-6; % [m]
N_CLAD = 1.45;
DELTA = 0.004;
N_CORE = N_CLAD * (1.0 + DELTA);
k0 = 2 * pi / WAVELENGTH;

% --- 2. 結合係数の設定 (論文の物理量に基づく) ---
% 隣接コアと対角コアで距離が異なることを考慮
pos = [ CORE_PITCH/2,  CORE_PITCH/2; -CORE_PITCH/2,  CORE_PITCH/2; ...
       -CORE_PITCH/2, -CORE_PITCH/2;  CORE_PITCH/2, -CORE_PITCH/2];
K_mn = zeros(NUM_CORES, NUM_CORES);
for m = 1:NUM_CORES
    for n = 1:NUM_CORES
        if m == n, continue; end
        dist = norm(pos(m,:) - pos(n,:));
        % 結合強度を調整 (論文のオーダーに合わせる)
        K_mn(m,n) = 1e-4 * exp(-dist / 8e-6); 
    end
end

% --- 3. シミュレーション (曲げ直径依存性) ---
bend_diameters_mm = logspace(2, 3.3, 50); 
results = zeros(length(bend_diameters_mm), 6); 
for i = 1:length(bend_diameters_mm)
    Rb = (bend_diameters_mm(i) / 2) * 1e-3; % 曲げ半径 [m]
    P = [1.0, 1e-9, 1e-9, 1e-9]; % コア1にのみ光を入力
    
    for z = 0:DZ:FIBER_LENGTH
        dP_dz = zeros(1, NUM_CORES);
        for m = 1:NUM_CORES
            for n = 1:NUM_CORES
                if m == n, continue; end
                % 式(12): 局所伝搬定数差
                delta_beta = k0 * (N_CORE - N_CLAD);
                delta_beta_prime = delta_beta + (CORE_PITCH / Rb);
                
                % 式(24): TAFによるパワー結合係数 (PCC)
                X = delta_beta_prime * DC_CORR / 2.0;
                h_mn = (K_mn(m,n)^2 * DC_CORR * sin(X+eps)^2) / (2.0 * X^2 + eps);
                
                dP_dz(m) = dP_dz(m) + h_mn * (P(n) - P(m));
            end
        end
        P = P + dP_dz * DZ;
    end
    
    % クロストーク計算 [dB]
    pairs = [1,2; 1,3; 1,4; 2,3; 2,4; 3,4];
    for p = 1:6
        results(i, p) = 10 * log10(max(P(pairs(p,2)), 1e-12) / P(pairs(p,1)));
    end
end

% --- 4. プロット ---
figure('Color', 'w', 'Position', [100, 100, 1000, 700]);
pair_names = {'1-2', '1-3', '1-4', '2-3', '2-4', '3-4'};
for p = 1:6
    subplot(2, 3, p);
    semilogx(bend_diameters_mm, results(:, p), 'LineWidth', 2);
    grid on;
    title(['Pair ', pair_names{p}]);
    xlabel('Bending diameter [mm]'); ylabel('Crosstalk [dB]');
    ylim([-80, -20]); xlim([100, 2000]);
end