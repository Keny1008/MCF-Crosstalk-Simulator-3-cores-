function acoplador_direcional_interativo()
    clc; close all;
    disp('===   SIMULATION OF TWO CORES ===');
    
    %% 1. パラメータ入力
    n_clad = input('クラッドの屈折率を入力: ');
    n_core1 = input('コアの屈折率 (n_1) を入力: ');
    n_core2 = n_core1;
    d1 = input('コア1の直径 [μm]を入力: ') * 1e-6;
    d2 = input('コア2の直径 [μm]を入力: ') * 1e-6;
    gap_choice = input('計算したいコア間隔 (Gap) [μm]: ') * 1e-6;
    e0 = 8.85418781e-12;             
    m0 = 1.25663706e-6;              
    c_speed = 1.0 / sqrt(e0 * m0);   
    f = 193e12;                      % 周波数 [Hz]
    lambda = c_speed / f;
    k0 = 2*pi/lambda;
    
    
    %% 2. 孤立モードの伝搬定数 (Beta) 計算
    f_iso1 = @(b) searching_beta(b, k0, n_clad, n_core1, d1);
    beta1 = solve_mode_robust(f_iso1, n_clad*k0, n_core1*k0);
    f_iso2 = @(b) searching_beta(b, k0, n_clad, n_core2, d2);
    beta2 = solve_mode_robust(f_iso2, n_clad*k0, n_core2*k0);
    
    Delta_beta = (beta1 - beta2) / 2;
    beta_avg = (beta1 + beta2) / 2;
    
    %% 3. 結合係数 (Kappa) の計算
    % 選択されたGapでの電界分布を取得してKappaを算出
    x_pert = linspace(-20e-6, 20e-6, 2000); % 広めの範囲で計算
    center1 = -gap_choice/2 - d1/2; 
    center2 = gap_choice/2 + d2/2; 
    
    E1 = calc_field_single_shifted(beta1, k0, n_clad, n_core1, d1, center1, x_pert);
    E2 = calc_field_single_shifted(beta2, k0, n_clad, n_core2, d2, center2, x_pert);
    
    I12 = (n_core2^2 - n_clad^2) * integral_overlap(x_pert, E1, E2, center2, d2);
    I21 = (n_core1^2 - n_clad^2) * integral_overlap(x_pert, E2, E1, center1, d1);
    
    Kappa = sqrt(((k0^2 / (2 * beta1)) * I12 / trapz(x_pert, E1.^2)) * ...
                 ((k0^2 / (2 * beta2)) * I21 / trapz(x_pert, E2.^2)));
             
    %% 4. パワービート (クロストーク解析) の計算とプロット
    Gamma = sqrt(Kappa^2 + Delta_beta^2);
    % P2(z) = (Kappa^2 / Gamma^2) * sin^2(Gamma * z)
    % クロストーク値 (デシベル): -10 * log10(P2_max)
    P2_max = (Kappa / Gamma)^2;
    crosstalk_dB = -10 * log10(P2_max);
    
    fprintf('\n--- 解析結果 ---\n');
    fprintf('結合係数 (Kappa): %.4e 1/m\n', Kappa);
    fprintf('最大クロストーク (P2_max): %.4f (%.2f dB)\n', P2_max, crosstalk_dB);
    
    % プロット
    % --- 修正後のプロット用コード ---
    % 距離を 0 から 10km (10,000m) まで広げる
    z = linspace(0, 1000, 100); 
    
    % クロストークの計算 (dB)
    % P2(z) = (Kappa^2 / Gamma^2) * sin^2(Gamma * z)
    % 10*log10(P2) でdB表示
    xtalk_dB = 10 * log10( (Kappa^2 / Gamma^2) * sin(Gamma * z).^2 );
    
    % -100dB以下の値はグラフ表示が乱れるので、クリップする
    xtalk_dB(xtalk_dB < -100) = -100;

    figure('Color', 'w');
    plot(z / 1000, xtalk_dB, 'b-', 'LineWidth', 2); % 単位をkmに変換
    title(sprintf('クロストーク特性 (Gap=%.1f μm)', gap_choice*1e6));
    xlabel('ファイバ長 [km]');
    ylabel('クロストーク [dB]');
    grid on;
    ylim([-100 0]); % 見やすいようにdBの範囲を調整    
end


function da_dz = ode_cmt_coupled_perturbation(~, a, beta_avg, Delta_beta, Kappa)
% Função EDO para CMT Perturbação (Assimétrica)
    
    da_dz = [
        -1i * (beta_avg + Delta_beta) * a(1) - 1i * Kappa * a(2);
        -1i * (beta_avg - Delta_beta) * a(2) - 1i * Kappa * a(1)
    ];
end
function I = integral_overlap(x_pert, E_i, E_j, center_j, d_j)
% Calcula a integral de superposição para Kappa: \int E_i * E_j * I(x) dx
    
    x_L = center_j - d_j / 2;
    x_R = center_j + d_j / 2;
    
    Perturbation_window = (x_pert >= x_L) & (x_pert <= x_R);
    
    Integrand = E_i .* E_j .* Perturbation_window;
    
    I = trapz(x_pert, Integrand);
end
function plot_field_profiles_comparative(k0, n_clad, n_core1, n_gap, n_core2, d1, d2, gaps_list, beta_min, beta_max)
% GRÁFICO 1 (Perfís Compostos vs. Isolados) - PLOTANDO CAMPO E_y (4 subplots)
    figure('Name', 'Fig 1: Perfis de Campo Ey (Assimetria) vs. Gap', 'Color', 'w', 'Position', [100 100 1000 800]);
    x_plot = linspace(-4e-6, 4e-6, 1000); 
    
    f_iso1 = @(b) error_func_single_slab(b, k0, n_clad, n_core1, d1);
    beta_iso1 = solve_mode_robust(f_iso1, beta_min, beta_max);
    f_iso2 = @(b) error_func_single_slab(b, k0, n_clad, n_core2, d2);
    beta_iso2 = solve_mode_robust(f_iso2, beta_min, beta_max);
    
    % Se algum modo fundamental não for encontrado, a função não plota.
    if isnan(beta_iso1) || isnan(beta_iso2)
        disp('AVISO: Modos Isolados não encontrados para plotagem da Fig 1.');
        return;
    end

    for i = 1:length(gaps_list)
        S_gap = gaps_list(i);
        subplot(2, 2, i); hold on;
        
        center_guide1 = -S_gap/2 - d1/2; 
        center_guide2 = S_gap/2 + d2/2;
        
        % Coordenadas das Interfaces
        x_g1_L = -S_gap/2 - d1; 
        x_g1_R = -S_gap/2; 
        x_g2_L = S_gap/2; 
        x_g2_R = S_gap/2 + d2;
        
        betas_sys = find_exact_modes_5layer(k0, n_clad, n_core1, n_gap, n_core2, d1, d2, S_gap, beta_min, beta_max);
        
        if length(betas_sys) < 2
            title(sprintf('Gap = %.1f \\mu m (Sem modos guiados)', S_gap*1e6)); 
            xlim([-3 3]); ylim([-0.5 1.2]); grid on;
            continue;
        end
        betas_sys = sort(betas_sys, 'descend'); 
        
        Ey_compound1 = calc_field_5layer_final(betas_sys(1), k0, n_clad, n_core1, n_gap, n_core2, d1, d2, S_gap, x_plot);
        Ey_compound2 = calc_field_5layer_final(betas_sys(2), k0, n_clad, n_core1, n_gap, n_core2, d1, d2, S_gap, x_plot);
        
        Ey_iso1 = calc_field_single_shifted(beta_iso1, k0, n_clad, n_core1, d1, center_guide1, x_plot);
        Ey_iso2 = calc_field_single_shifted(beta_iso2, k0, n_clad, n_core2, d2, center_guide2, x_plot);
        
        % Lógica de Associação e Normalização
        overlap1_vs_iso1 = trapz(x_plot, Ey_compound1 .* Ey_iso1); 
        overlap1_vs_iso2 = trapz(x_plot, Ey_compound1 .* Ey_iso2); 
        Comp1_is_Guia1 = abs(overlap1_vs_iso1) > abs(overlap1_vs_iso2);
        
        if Comp1_is_Guia1
             Ey_Guia1_Comp = Ey_compound1; Ey_Guia2_Comp = Ey_compound2;
        else
             Ey_Guia1_Comp = Ey_compound2; Ey_Guia2_Comp = Ey_compound1; 
        end
        
        [~, idx_c1] = min(abs(x_plot - center_guide1));
        [~, idx_c2] = min(abs(x_plot - center_guide2));
        
        Ey_iso1_norm = Ey_iso1 / Ey_iso1(idx_c1); 
        Ey_Guia1_Comp_norm = Ey_Guia1_Comp / Ey_Guia1_Comp(idx_c1);
        if sign(Ey_iso1_norm(idx_c1)) ~= sign(Ey_Guia1_Comp_norm(idx_c1)), Ey_Guia1_Comp_norm = -Ey_Guia1_Comp_norm; end
        
        Ey_iso2_norm = Ey_iso2 / Ey_iso2(idx_c2); 
        Ey_Guia2_Comp_norm = Ey_Guia2_Comp / Ey_Guia2_Comp(idx_c2);
        if sign(Ey_iso2_norm(idx_c2)) ~= sign(Ey_Guia2_Comp_norm(idx_c2)), Ey_Guia2_Comp_norm = -Ey_Guia2_Comp_norm; end
        
        plot(x_plot*1e6, Ey_Guia1_Comp_norm, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Comp. (Guia 1)'); 
        plot(x_plot*1e6, Ey_iso1_norm, 'k:', 'LineWidth', 2.0, 'DisplayName', 'Isolado 1');
        plot(x_plot*1e6, Ey_Guia2_Comp_norm, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Comp. (Guia 2)'); 
        plot(x_plot*1e6, Ey_iso2_norm, 'k--', 'LineWidth', 2.0, 'DisplayName', 'Isolado 2');       
        
        % Linhas Verticais com Legendas de Interfaces
        xline(x_g1_L*1e6, 'k:', 'LineWidth', 0.5, 'HandleVisibility', 'off'); 
        xline(x_g1_R*1e6, 'k:', 'LineWidth', 0.5, 'DisplayName', 'Int. Guia 1 / Gap');
        xline(x_g2_L*1e6, 'k:', 'LineWidth', 0.5, 'DisplayName', 'Int. Gap / Guia 2');
        xline(x_g2_R*1e6, 'k:', 'LineWidth', 0.5, 'HandleVisibility', 'off'); 
        
        title(sprintf('Gap = %.1f \\mu m (d1=%.4f, d2=%.4f \\mu m)', S_gap*1e6, d1*1e6, d2*1e6));
        xlim([-3 3]); ylim([-0.5 1.2]); grid on;
        
        if i == 1
            legend('show', 'Location', 'northeast');
            ylabel('E_y Field');
        end
        if i > 2
            xlabel('x (\\mu m)');
        end
    end
end
function [Ey] = calc_field_single_shifted(beta, k0, n_clad, n_core, d, center_x, x_vec)
    % Campo do guia isolado transladado
    gamma = sqrt(beta^2 - n_clad^2*k0^2);
    kappa = sqrt(n_core^2*k0^2 - beta^2);
    Ey = zeros(size(x_vec));
    
    for k = 1:length(x_vec)
        x_local = x_vec(k) - center_x;
        
        if abs(x_local) <= d/2
            Ey(k) = cos(kappa * x_local);
        else
            dist = abs(x_local) - d/2;
            Ey(k) = cos(kappa * d/2) * exp(-gamma * dist);
        end
    end
end
function [Ey] = calc_field_5layer_final(beta, k0, n1, n2, n3, n4, d1, d2, gap, x_vec)
    % Campo exato da estrutura de 5 camadas (TMM)
    gamma1 = sqrt(beta^2 - n1^2*k0^2); kappa2 = sqrt(n2^2*k0^2 - beta^2);
    gamma3 = sqrt(beta^2 - n3^2*k0^2); kappa4 = sqrt(n4^2*k0^2 - beta^2);
    
    if ~isreal(gamma1) || ~isreal(kappa2) || ~isreal(kappa4), Ey = NaN*x_vec; return; end
    
    Ey = zeros(size(x_vec));
    x1 = -gap/2 - d1; x2 = -gap/2; x3 = gap/2; x4 = gap/2 + d2;
    A = 1.0; 
    
    for i = 1:length(x_vec)
        x = x_vec(i);
        
        if x < x1
            Ey(i) = A * exp(gamma1 * (x - x1));
        
        elseif x >= x1 && x < x2
            dx = x - x1;
            Ey(i) = A * cos(kappa2*dx) + (A*gamma1/kappa2) * sin(kappa2*dx);
        
        elseif x >= x2 && x < x3
            E_x2 = A*cos(kappa2*d1) + (A*gamma1/kappa2)*sin(kappa2*d1);
            dE_x2 = -A*kappa2*sin(kappa2*d1) + A*gamma1*cos(kappa2*d1);
            dx = x - x2;
            Ey(i) = E_x2 * cosh(gamma3*dx) + (dE_x2/gamma3) * sinh(gamma3*dx);
        
        elseif x >= x3 && x < x4
            E_x2 = A*cos(kappa2*d1) + (A*gamma1/kappa2)*sin(kappa2*d1); dE_x2 = -A*kappa2*sin(kappa2*d1) + A*gamma1*cos(kappa2*d1);
            E_x3 = E_x2*cosh(gamma3*gap) + (dE_x2/gamma3)*sinh(gamma3*gap); dE_x3 = E_x2*gamma3*sinh(gamma3*gap) + dE_x2*cosh(gamma3*gap);
            dx = x - x3;
            Ey(i) = E_x3 * cos(kappa4*dx) + (dE_x3/kappa4) * sin(kappa4*dx);
        
        else
            E_x2 = A*cos(kappa2*d1) + (A*gamma1/kappa2)*sin(kappa2*d1); dE_x2 = -A*kappa2*sin(kappa2*d1) + A*gamma1*cos(kappa2*d1);
            E_x3 = E_x2*cosh(gamma3*gap) + (dE_x2/gamma3)*sinh(gamma3*gap); dE_x3 = E_x2*gamma3*sinh(gamma3*gap) + dE_x2*cosh(gamma3*gap);
            E_x4 = E_x3*cos(kappa4*d2) + (dE_x3/kappa4)*sin(kappa4*d2);
            Ey(i) = E_x4 * exp(-gamma1 * (x - x4));
        end
    end
end
function betas = find_exact_modes_5layer(k0, n1, n2, n3, n4, d1, d2, gap, b_min, b_max)
    % Solver TMT exato para 5 camadas
    f_res = @(b) solve_boundary_error(b, k0, n1, n2, n3, n4, d1, d2, gap);
    b_step = linspace(b_min, b_max, 500);
    res_vals = arrayfun(f_res, b_step);
    betas = [];
    for k = 1:length(b_step)-1
        if sign(res_vals(k)) ~= sign(res_vals(k+1))
            try
                root = fzero(f_res, [b_step(k), b_step(k+1)]);
                betas = [betas, root];
            catch; end
        end
    end
end
function residual = solve_boundary_error(beta, k0, n1, n2, n3, n4, d1, d2, gap)
    % Erro de contorno para 5 camadas (TE)
    gamma1 = sqrt(beta^2 - n1^2*k0^2); kappa2 = sqrt(n2^2*k0^2 - beta^2);
    gamma3 = sqrt(beta^2 - n3^2*k0^2); kappa4 = sqrt(n4^2*k0^2 - beta^2);
    
    if ~isreal(gamma1) || ~isreal(kappa2) || ~isreal(kappa4), residual = NaN; return; end
    
    M_core1 = [cos(kappa2*d1), sin(kappa2*d1)/kappa2; -kappa2*sin(kappa2*d1), cos(kappa2*d1)];
    M_gap   = [cosh(gamma3*gap), sinh(gamma3*gap)/gamma3; gamma3*sinh(gamma3*gap), cosh(gamma3*gap)];
    M_core2 = [cos(kappa4*d2), sin(kappa4*d2)/kappa4; -kappa4*sin(kappa4*d2), cos(kappa4*d2)];
    
    State = [1; gamma1];
    State = M_core2 * M_gap * M_core1 * State; 
    
    E_final = State(1); dE_final = State(2);
    residual = dE_final + gamma1 * E_final; 
end


function res = searching_beta(beta, k0, n_clad, n_core, d)
    % Equação Característica do Slab Único (Modo TE Simétrico)
    gamma = sqrt(beta^2 - n_clad^2*k0^2);
    kappa = sqrt(n_core^2*k0^2 - beta^2);
    res = kappa * d - 2*atan(gamma/kappa); 
end


function beta_sol = solve_mode_robust(func_handle, b_min, b_max)
    % Solver robusto para encontrar betas
    b_vec = linspace(b_min, b_max, 200);
    res_vec = arrayfun(func_handle, b_vec);
    beta_sol = NaN;
    for i = length(b_vec)-1:-1:1
        if sign(res_vec(i)) ~= sign(res_vec(i+1))
            try
                beta_sol = fzero(func_handle, [b_vec(i), b_vec(i+1)]);
                return;
            catch; end
        end
    end
end