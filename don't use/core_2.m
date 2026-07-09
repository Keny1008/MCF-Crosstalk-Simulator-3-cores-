function acoplador_direcional_interativo()
    clc; close all;
    disp('===SIMULATION OF TWO CORES===');
    
    %% 1. パラメータ
    n_clad = input('クラッドの屈折率を入力: ');
    n_core1 = input('コア1の屈折率 (n_1) を入力: ');
    n_core2 = input('コア2の屈折率 (n_2) を入力: ');
    d1 = input('コア1の直径 [μm]を入力: ') * 1e-6;
    d2 = input('コア2の直径 [μm]を入力: ') * 1e-6;
    gap_choice = input('結合解析用のギャップ間隔 [μm] (S_um): ') * 1e-6;

    e0 = 8.85418781e-12;             
    m0 = 1.25663706e-6;              
    c_speed = 1.0 / sqrt(e0 * m0);   
    f = 193e12;                      % 周波数 [Hz]
    lambda = c_speed / f;
    k0 = 2*pi/lambda;
    
    % Lista de Gaps para análise de perfis e Kappa
    gaps_list = [0.2e-6, 0.5e-6, 1.0e-6, 2.0e-6]; 
    if ~ismember(gap_choice, gaps_list)
        gaps_list = sort(unique([gaps_list, gap_choice]));
    end
    
    n_gap = n_clad; % O Gap é o mesmo material do Cladding
    
    % Range de busca do Beta (TE0 - Modo Fundamental)
    beta_min = n_clad * k0 + 100; 
    beta_max = max([n_core1, n_core2]) * k0 - 100;
    
    %% 2. CÁLCULO DE MODOS ISOLADOS E PARÂMETROS CMT
    fprintf('\n--- 2. Cálculo de Modos Isolados ---\n');
    
    f_iso1 = @(b) error_func_single_slab(b, k0, n_clad, n_core1, d1);
    beta1 = solve_mode_robust(f_iso1, beta_min, beta_max);
    
    f_iso2 = @(b) error_func_single_slab(b, k0, n_clad, n_core2, d2);
    beta2 = solve_mode_robust(f_iso2, beta_min, beta_max);
    
    if isnan(beta1) || isnan(beta2)
         error('ERRO: Não foi possível encontrar os modos fundamentais (TE0) isolados. Ajuste os parâmetros de guia.');
    end
    
    fprintf('  Beta Isolado 1: %.4e 1/m\n', beta1);
    fprintf('  Beta Isolado 2: %.4e 1/m\n', beta2);
    
    Delta_beta = (beta1 - beta2) / 2;
    beta_avg_pert = (beta1 + beta2) / 2;
    
    fprintf('  Dessintonia (Delta_beta): %.4e 1/m\n', Delta_beta);
    
    
    %% 3. GRÁFICO 1: Perfis de Campo Ey (Comparativo da Assimetria)
    plot_field_profiles_comparative(k0, n_clad, n_core1, n_gap, n_core2, d1, d2, gaps_list, beta_min, beta_max);
    
    
    %% 4. GRÁFICO 2: Coeficiente de Acoplamento (Kappa vs. Gap)
    fprintf('\n--- 4. Cálculo de Kappa vs. Gap (CMT Perturbação) ---\n');
    
    Kappa_list = zeros(size(gaps_list));
    x_pert = linspace(-4.0e-6, 4.0e-6, 2000);
    
    for i = 1:length(gaps_list)
        S_gap = gaps_list(i);
        
        center1 = -S_gap/2 - d1/2; center2 = S_gap/2 + d2/2; 
        
        E1_campo = calc_field_single_shifted(beta1, k0, n_clad, n_core1, d1, center1, x_pert);
        E2_campo = calc_field_single_shifted(beta2, k0, n_clad, n_core2, d2, center2, x_pert);
        
        V11_iso = trapz(x_pert, E1_campo.^2);
        V22_iso = trapz(x_pert, E2_campo.^2);
        
        I_12 = (n_core2^2 - n_clad^2) * integral_overlap(x_pert, E1_campo, E2_campo, center2, d2);
        I_21 = (n_core1^2 - n_clad^2) * integral_overlap(x_pert, E2_campo, E1_campo, center1, d1);
        
        K12 = (k0^2 / (2 * beta1)) * I_12 / V11_iso;
        K21 = (k0^2 / (2 * beta2)) * I_21 / V22_iso;
        
        Kappa_list(i) = sqrt(K12 * K21); 
        
        fprintf('  Gap=%.1f um: Kappa=%.4e 1/m\n', S_gap*1e6, Kappa_list(i));
    end
    
    figure('Name', 'Fig 2: Coeficiente de Acoplamento (Kappa) vs. Gap', 'Color', 'w');
    plot(gaps_list * 1e6, Kappa_list, 'k-o', 'LineWidth', 2, 'MarkerSize', 8);
    title('Coeficiente de Acoplamento \kappa (CMT Perturbação) vs. Separação');
    xlabel('Gap S (\mu m)');
    ylabel('\kappa (m^{-1})');
    grid on; 
    
    
    %% 5. GRÁFICO 3: Batimento de Potência (EDOs da CMT)
    
    Gap_edo = gap_choice; 
    
    [~, idx_gap_edo] = min(abs(gaps_list - Gap_edo));
    Kappa_edo = Kappa_list(idx_gap_edo);

    Gamma_pert = sqrt(Kappa_edo^2 + Delta_beta^2);
    
    fprintf('\n--- 5. Solução EDO (Gap = %.4f um) ---\n', Gap_edo*1e6);
    
    Lc_pert_m = pi / (2 * Gamma_pert);
    Lc_pert_um = Lc_pert_m * 1e6;
    
    f_ode = @(z, a) ode_cmt_coupled_perturbation(z, a, beta_avg_pert, Delta_beta, Kappa_edo);
    A0 = [1; 0]; 
    
    Z_span = [0, 4 * Lc_pert_m];
    
    disp('Resolvendo EDOs da CMT Perturbação...');
    [Z, A] = ode45(f_ode, Z_span, A0); 
    
    P1 = abs(A(:, 1)).^2;
    P2 = abs(A(:, 2)).^2;
    
    % --- GRÁFICO 3: Batimento de Potência ---
    figure('Name', 'Fig 3: Batimento de Potência (CMT)', 'Color', 'w');
    plot(Z * 1e6, P1, 'b-', 'LineWidth', 2, 'DisplayName', 'P1 (Guia 1)');
    hold on;
    plot(Z * 1e6, P2, 'r--', 'LineWidth', 2, 'DisplayName', 'P2 (Guia 2)');
    
    P_max_transfer = (Kappa_edo / Gamma_pert)^2; 
    
    title(sprintf('Batimento de Potência CMT-Perturb.: Gap=%.4f \\mu m (Lc=%.2f \\mu m)', Gap_edo*1e6, Lc_pert_um));
    xlabel('Comprimento de Interação z (\\mu m)');
    ylabel('Potência Normalizada |a_j|^2');
    yline(P_max_transfer, 'k:', 'LineWidth', 1.0, 'DisplayName', 'P_{max, transferido}');
    legend('show', 'Location', 'northeast');
    grid on; ylim([0, 1.05]); 
    
    fprintf('  Comprimento de Acoplamento (Lc): %.2f um\n', Lc_pert_um);
    fprintf('  Potência Máx. Transferida (P2max): %.3f\n', P_max_transfer);
    
    
    %% 6. GRÁFICO 4: Modos Compostos (TMM - Ey) no Gap Escolhido
    
    Gap_fig = gap_choice; 
    center_guide1 = -Gap_fig/2 - d1/2; 
    center_guide2 = Gap_fig/2 + d2/2;

    % Coordenadas das Interfaces para Legenda
    x_g1_L = -Gap_fig/2 - d1; 
    x_g1_R = -Gap_fig/2; % Interface Guia 1 / Gap
    x_g2_L = Gap_fig/2; % Interface Gap / Guia 2
    x_g2_R = Gap_fig/2 + d2;

    betas_sys_fig = find_exact_modes_5layer(k0, n_clad, n_core1, n_gap, n_core2, d1, d2, Gap_fig, beta_min, beta_max);
    
    fprintf('\n--- 6. Análise TMM (Gap = %.4f um) ---\n', Gap_fig*1e6);
    
    if length(betas_sys_fig) < 2
        fprintf('AVISO: Não foi possível plotar Fig. 4. Menos de 2 modos guiados (TE0) encontrados para Gap=%.4f um.\n', Gap_fig*1e6);
    else
        betas_sys_fig = sort(betas_sys_fig, 'descend'); 
        
        beta_compound1 = betas_sys_fig(1); % Modo Par (Maior Beta)
        beta_compound2 = betas_sys_fig(2); % Modo Ímpar (Segundo Maior Beta)
        
        fprintf('  Beta Modo Par (TMM): %.4e 1/m\n', beta_compound1);
        fprintf('  Beta Modo Ímpar (TMM): %.4e 1/m\n', beta_compound2);
        
        x_plot_long = linspace(-4e-6, 4e-6, 2000); 
        Ey_compound1 = calc_field_5layer_final(beta_compound1, k0, n_clad, n_core1, n_gap, n_core2, d1, d2, Gap_fig, x_plot_long);
        Ey_compound2 = calc_field_5layer_final(beta_compound2, k0, n_clad, n_core1, n_gap, n_core2, d1, d2, Gap_fig, x_plot_long);

        % Normalização e Ajuste de Fase
        Max_Ey = max(max(abs(Ey_compound1)), max(abs(Ey_compound2)));
        Ey_norm1 = Ey_compound1 / Max_Ey;
        Ey_norm2 = Ey_compound2 / Max_Ey;

        % Ajuste de sinal para garantir Modo Par e Modo Ímpar visualmente
        [~, idx_c1] = min(abs(x_plot_long - center_guide1));
        if Ey_norm1(idx_c1) < 0, Ey_norm1 = -Ey_norm1; end
        if Ey_norm2(idx_c1) > 0, Ey_norm2 = -Ey_norm2; end 

        figure('Name', 'Fig 4: Modos Compostos (Ey) no Gap Escolhido (TMM)', 'Color', 'w');
        plot(x_plot_long * 1e6, Ey_norm1, 'k-', 'LineWidth', 2.0, 'DisplayName', 'Modo Par (\beta_e)');
        hold on;
        plot(x_plot_long * 1e6, Ey_norm2, 'k:', 'LineWidth', 2.0, 'DisplayName', 'Modo Ímpar (\beta_o)'); 
        
        % Linhas Verticais com Legendas de Interfaces
        xline(x_g1_R*1e6, 'k:', 'LineWidth', 0.5, 'DisplayName', 'Int. Guia 1 / Gap');
        xline(x_g2_L*1e6, 'k:', 'LineWidth', 0.5, 'DisplayName', 'Int. Gap / Guia 2');
        
        xline(x_g1_L*1e6, 'k:', 'LineWidth', 0.5, 'HandleVisibility', 'off'); 
        xline(x_g2_R*1e6, 'k:', 'LineWidth', 0.5, 'HandleVisibility', 'off'); 

        title(sprintf('Modos Compostos TE (Ey): Gap = %.4f \\mu m', Gap_fig*1e6));
        xlabel('X-Axis (\\mu m)');
        ylabel('E_y(x) (Normalizado)');
        grid on; xlim([-2.5 2.5]); ylim([-1.1 1.1]);
        legend('show', 'Location', 'northeast');
    end
end

% =========================================================================
% === FUNÇÕES AUXILIARES (Devem ser mantidas como funções locais no MATLAB) ===
% =========================================================================

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
function res = error_func_single_slab(beta, k0, n_clad, n_core, d)
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