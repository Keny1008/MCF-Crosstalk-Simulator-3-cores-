function drag_interactive_mcf_crosstalk_Ncore()
    % =========================================================================
    % マルチコアファイバ,クロストーク解析 (任意コア数 N 対応版・クラッド固定)
    % =========================================================================
    
    % --- 1. 固定の物理パラメータ設定 ---
    e0 = 8.85418781e-12;             % 真空の誘電率
    m0 = 1.25663706e-6;              % 真空の透磁率
    c_speed = 1.0 / sqrt(e0 * m0);   % 光速
    f = 193e12;                      % 周波数 [Hz]
    lambda = c_speed / f;            % 波長    
    n_clad = input('クラッドの屈折率を入力: ');
    n_core = input('コアの屈折率を入力: ');
    a = input('コアの半径 [μm]を入力: ') * 1e-6;
    N = input('コア数を入力: ');     % コア数を可変に
    FIBER_LENGTH = 10000;          
    SEGMENT_LEN = 1.0;             
    NUM_SEGMENTS = FIBER_LENGTH / SEGMENT_LEN;
    z = linspace(0, FIBER_LENGTH, NUM_SEGMENTS+1);

    % --- 2. GUIの作成 ---
    fig = uifigure('Name', 'Independent Draggable MCF Crosstalk (N-core)', 'Position', [100, 100, 1150, 550]);
    gl = uigridlayout(fig, [2, 2]);
    gl.RowHeight = {'1x', 40};
    gl.ColumnWidth = {'1.5x', '1x'};
    
    % クロストークグラフ用Axes
    ax_graph = uiaxes(gl);
    ax_graph.Layout.Row = 1; ax_graph.Layout.Column = 1;
    xlabel(ax_graph, 'Distance [m]'); ylabel(ax_graph, 'Crosstalk [dB]');
    grid(ax_graph, 'on'); hold(ax_graph, 'on');
    
    % 断面図用Axes
    ax_cross = uiaxes(gl);
    ax_cross.Layout.Row = 1; ax_cross.Layout.Column = 2;
    xlabel(ax_cross, 'x [m]'); ylabel(ax_cross, 'y [m]');
    axis(ax_cross, 'equal'); grid(ax_cross, 'on'); hold(ax_cross, 'on');
    title(ax_cross, 'すべてのコアを自由にドラッグして配置できます');
    
    % ※表示範囲をクラッドが綺麗に収まるように対称に固定
    xlim(ax_cross, [-100e-6, 100e-6]); ylim(ax_cross, [-100e-6, 100e-6]); 
    
    % 情報表示用ラベル
    lbl_info = uilabel(gl, 'Text', '初期化中...', 'FontWeight', 'bold', 'FontSize', 12, ...
                        'WordWrap', 'on');
    lbl_info.Layout.Row = 2; lbl_info.Layout.Column = [1 2];

    % --- 3. モード解析と結合係数の事前計算 ---
    wb = uiprogressdlg(fig, 'Title', '初期化中', 'Message', '各コアの距離変動に対応するテーブルを作成しています...');
    
    NA = sqrt(n_core^2 - n_clad^2);
    V_param = (2 * pi * a / lambda) * NA;
    char_eq = @(u) (besselj(0, u) ./ (u .* besselj(1, u))) - ...
                   (besselk(0, sqrt(max(eps, V_param^2 - u^2))) ./ (sqrt(max(eps, V_param^2 - u^2)) .* besselk(1, sqrt(max(eps, V_param^2 - u^2)))));
    u_val = fzero(char_eq, V_param/2); 
    w_val = sqrt(max(eps, V_param^2 - u_val^2));
    beta = sqrt((2 * pi * n_core / lambda)^2 - (u_val/a)^2);
    
    N_norm = 1 / sqrt(integral(@(r) 2*pi*r .* ((r < a) .* (besselj(0, u_val*r/a)/besselj(0, u_val)).^2 + (r >= a) .* (besselk(0, w_val*r/a)/besselk(0, w_val)).^2), 0, 10*a));
    E_norm = @(r) N_norm * ((r < a) .* (besselj(0, u_val*r/a)/besselj(0, u_val)) + (r >= a) .* (besselk(0, w_val*r/a)/besselk(0, w_val)));
    
    % 距離が変化しても対応できるよう、11um ~ 150um までの広範囲でkappaを事前計算
    pitch_array = linspace(11e-6, 150e-6, 60);
    kappa_array = zeros(size(pitch_array));
    for i = 1:length(pitch_array)
        p_temp = pitch_array(i);
        integrand = @(x, y) (n_core^2 - n_clad^2) .* E_norm(sqrt(x.^2 + y.^2)) .* E_norm(sqrt((x-p_temp).^2 + y.^2));
        kappa_array(i) = real( ( (2 * pi * f * e0) / (2 * beta) ) * integral2(integrand, p_temp-a, p_temp+a, -a, a) );
        wb.Value = i / length(pitch_array);
    end
    close(wb);

    % --- 4. コア座標の初期設定 (N個をリング状に自動配置) ---
    p_init = 39.2e-6;
    core_pos = zeros(N, 2);
    core_pos(1,:) = [0, 0]; % Core 1 を原点に
    if N > 1
        theta = linspace(0, 2*pi, N)';   
        theta = theta(1:end-1);          
        core_pos(2:end, 1) = p_init * cos(theta);
        core_pos(2:end, 2) = p_init * sin(theta);
    end

    % --- 5. 描画オブジェクトの初期化 ---
    % 【修正点】クラッドを中心(0,0)に固定し、半径も固定して描画
    fixed_clad_r = p_init + 5*a; % 余裕をもたせたクラッド半径 (約64μm)
    h_clad = rectangle(ax_cross, 'Position', [-fixed_clad_r, -fixed_clad_r, 2*fixed_clad_r, 2*fixed_clad_r], ...
                        'Curvature', [1 1], 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'k', 'PickableParts', 'none');
    
    colors = lines(N);           
    h_core  = gobjects(N,1);
    h_text  = gobjects(N,1);
    
    for k = 1:N
        col = colors(k,:);
        if k == 1
            col = [1 0 0]; % Core1(入力コア)は赤で強調
        end
        h_core(k) = rectangle(ax_cross, 'Position', [0 0 1 1], 'Curvature', [1 1], ...
                               'FaceColor', col, 'EdgeColor', 'k', 'PickableParts', 'all');
        h_text(k) = text(ax_cross, 0, 0, sprintf('Core %d', k), ...
                          'HorizontalAlignment', 'center', 'PickableParts', 'none', 'FontWeight', 'bold');
    end

    % グラフ用オブジェクト(N本)
    h_line = gobjects(N,1);
    leg_labels = cell(N,1);
    for k = 1:N
        col = colors(k,:);
        if k == 1
            col = [0 0 0]; % Core1は黒(Input)
            lw = 2;
        else
            lw = 1.5;
        end
        h_line(k) = plot(ax_graph, z, zeros(size(z)), 'Color', col, 'LineWidth', lw);
        if k == 1
            leg_labels{k} = 'Core 1 (Input)';
        else
            leg_labels{k} = sprintf('Core %d', k);
        end
    end
    ylim(ax_graph, [-160, 5]); xlim(ax_graph, [0, FIBER_LENGTH]);
    legend(ax_graph, h_line, leg_labels, 'Location', 'southwest');

    % --- 6. ドラッグ操作のコールバック設定 ---
    current_drag_core = 0; 
    
    for k = 1:N
        h_core(k).ButtonDownFcn = @(~,~) startDrag(k);
    end
    fig.WindowButtonMotionFcn = @dragging;
    fig.WindowButtonUpFcn = @stopDrag;

    % 初回描画
    updateAll();

    % --- 7. コールバック・更新関数 ---
    function startDrag(core_idx)
        current_drag_core = core_idx;
    end

    function dragging(~, ~)
        if current_drag_core == 0, return; end
        
        cp = ax_cross.CurrentPoint;
        core_pos(current_drag_core, 1) = cp(1,1);
        core_pos(current_drag_core, 2) = cp(1,2);
        
        updateAll(); 
    end

    function stopDrag(~, ~)
        current_drag_core = 0;
    end

    function updateAll()
        % 1. 全コアペア間の距離(ピッチ)行列を計算
        D = zeros(N, N);
        for ii = 1:N
            for jj = ii+1:N
                d_ij = norm(core_pos(ii,:) - core_pos(jj,:));
                D(ii,jj) = d_ij;
                D(jj,ii) = d_ij;
            end
        end
        
        % 2. 断面図の更新
        % 【修正点】クラッドのPosition更新処理を削除。コアとテキストのみ動かす。
        for k2 = 1:N
            set(h_core(k2), 'Position', [core_pos(k2,1) - a, core_pos(k2,2) - a, 2*a, 2*a]);
            set(h_text(k2), 'Position', [core_pos(k2,1), core_pos(k2,2) - 1.5*a]);
        end
        
        % 3. クロストークの瞬時計算
        p_min = min(pitch_array); p_max = max(pitch_array);
        
        M = zeros(N, N);
        for ii = 1:N
            for jj = ii+1:N
                d_clip = min(max(D(ii,jj), p_min), p_max);
                k_ij = interp1(pitch_array, kappa_array, d_clip, 'spline');
                M(ii,jj) = k_ij;
                M(jj,ii) = k_ij;
            end
        end
        
        [V_mat, D_mat] = eig(M);
        diag_D = diag(D_mat);
        A0 = zeros(N,1); A0(1) = 1;   
        V_inv_A0 = V_mat \ A0;
        
        kz_mat = -1j * diag_D * z; 
        A_all = V_mat * (diag(V_inv_A0) * exp(kz_mat)); 
        P = abs(A_all).^2; 
        
        % 4. グラフの更新
        for k2 = 1:N
            set(h_line(k2), 'YData', 10*log10(P(k2, :) + 1e-20));
        end
        
        % 5. 情報表示テキスト
        info_str = '';
        for ii = 1:N
            for jj = ii+1:N
                info_str = [info_str, sprintf('[C%d-C%d]: %.2f \\mum   ', ii, jj, D(ii,jj)*1e6)]; %#ok<AGROW>
            end
        end
        lbl_info.Text = ['ピッチ  ', info_str];
    end
end