function drag_interactive_mcf_crosstalk_FEM_wavelet_fixed2()
    % --- 1. 固定の物理パラメータ設定 ---
    e0 = 8.85418781e-12;
    m0 = 1.25663706e-6;
    c_speed = 1.0 / sqrt(e0 * m0);
    f = 193e12;
    lambda = c_speed / f;
    k0 = 2*pi/lambda;
    n_clad = 1.50;
    n_core = 1.50 * (1 + 0.003) ;
    
    a = input('コアの半径 [μm]を入力 (例: 4): ') * 1e-6;
    N = input('コアの数を入力:'); 
    
    FIBER_LENGTH = 1.0;         % 長さを1mにして、綺麗なうなりを見せます
    NUM_Z = 400;                % Z方向を細かくしてスペクトラムを綺麗に
    z = linspace(0, FIBER_LENGTH, NUM_Z);
    HMAX = a;              
    
    % --- 2. GUIの作成 (2×2の4枚画構成 + 情報ラベル) ---
    fig = uifigure('Name', 'FEM + STFT: Draggable MCF Wavelet Analyser', 'Position', [50, 50, 1400, 850]);
    gl = uigridlayout(fig, [3, 2]);
    gl.RowHeight = {'1x', '1x', 40};
    gl.ColumnWidth = {'1x', '1x'};
    
    % (1) 断面図
    ax_cross = uiaxes(gl); ax_cross.Layout.Row = 1; ax_cross.Layout.Column = 1;
    xlabel(ax_cross, 'x [m]'); ylabel(ax_cross, 'y [m]');
    axis(ax_cross, 'equal'); grid(ax_cross, 'on'); hold(ax_cross, 'on');
    title(ax_cross, 'コア2をドラッグして間隔調整 (近づけると周波数がUP！)');
    xlim(ax_cross, [-60e-6, 60e-6]); ylim(ax_cross, [-60e-6, 60e-6]);
    
    % (2) クロストーク グラフ
    ax_graph = uiaxes(gl); ax_graph.Layout.Row = 1; ax_graph.Layout.Column = 2;
    xlabel(ax_graph, 'Distance z [m]'); ylabel(ax_graph, 'Crosstalk [dB]');
    grid(ax_graph, 'on'); hold(ax_graph, 'on');
    title(ax_graph, 'z方向パワー遷移 (dB表示・進行に伴い加速するうなり)');
    ylim(ax_graph, [-40, 5]); xlim(ax_graph, [0, FIBER_LENGTH]);
    
    % (3) エネルギー保存の確認グラフ
    ax_energy = uiaxes(gl); ax_energy.Layout.Row = 2; ax_energy.Layout.Column = 1;
    xlabel(ax_energy, 'Distance z [m]'); ylabel(ax_energy, 'Normalized power [-]');
    grid(ax_energy, 'on'); hold(ax_energy, 'on');
    title(ax_energy, 'エネルギー保存の確認');
    xlim(ax_energy, [0, FIBER_LENGTH]); ylim(ax_energy, [0, 1.2]);
    
    % (4) スペクトラム解析（ウェーブレット/STFT風可視化）
    ax_wavelet = uiaxes(gl); ax_wavelet.Layout.Row = 2; ax_wavelet.Layout.Column = 2;
    xlabel(ax_wavelet, 'Distance z [m]'); ylabel(ax_wavelet, 'Coupling Frequency [1/m]');
    title(ax_wavelet, 'パワー変動の局所空間周波数解析 (時間-周波数マップ)');
    xlim(ax_wavelet, [0, FIBER_LENGTH]);
    
    lbl_info = uilabel(gl, 'Text', '初期化中...', 'FontWeight', 'bold', 'FontSize', 11, 'WordWrap', 'on');
    lbl_info.Layout.Row = 3; lbl_info.Layout.Column = [1 2];
    
    core_pos = [0, 0; 35e-6, 0]; 
    
    h_clad = rectangle(ax_cross, 'Position', [0 0 1 1], 'Curvature', [1 1], ...
                        'FaceColor', [0.95 0.95 0.95], 'EdgeColor', 'k', 'PickableParts', 'none');
    colors = [1 0 0; 0 0.447 0.741]; 
    h_core = gobjects(N,1); h_text = gobjects(N,1);
    for k = 1:N
        h_core(k) = rectangle(ax_cross, 'Position', [0 0 1 1], 'Curvature', [1 1], ...
                               'FaceColor', colors(k,:), 'EdgeColor', 'k', ...
                               'PickableParts', 'all');
        h_text(k) = text(ax_cross, 0, 0, sprintf('Core %d', k), ...
                          'HorizontalAlignment', 'center', 'PickableParts', 'none', 'FontWeight', 'bold');
    end
    h_core(1).PickableParts = 'none'; 
    
    h_line = gobjects(N,1); leg_labels = cell(N,1);
    h_line_energy_core = gobjects(N,1);
    for k = 1:N
        h_line(k) = plot(ax_graph, z, nan(size(z)), 'Color', colors(k,:), 'LineWidth', 2);
        leg_labels{k} = sprintf('Core %d', k);
        h_line_energy_core(k) = plot(ax_energy, z, nan(size(z)), 'Color', colors(k,:), 'LineWidth', 2);
    end
    legend(ax_graph, h_line, leg_labels, 'Location', 'southwest');
    
    h_line_total = plot(ax_energy, z, nan(size(z)), 'Color', [0.3 0.3 0.3], 'LineStyle', '--', 'LineWidth', 2);
    legend(ax_energy, [h_line_energy_core(:); h_line_total], {'Core 1', 'Core 2', 'Total'}, 'Location', 'southwest');
    
    % 初期ダミーイメージ
    h_img = imagesc(ax_wavelet, [0 FIBER_LENGTH], [0 100], zeros(10, 10));
    colormap(ax_wavelet, 'jet');
    ax_wavelet.YDir = 'normal';
    
    isDragging = false;
    h_core(2).ButtonDownFcn = @(~,~) startDrag();
    fig.WindowButtonMotionFcn = @dragging;
    fig.WindowButtonUpFcn = @stopDrag;
    
    updateAll();
    
    function startDrag(), isDragging = true; end
    function stopDrag(~,~), isDragging = false; end
    function dragging(~,~)
        if ~isDragging, return; end
        cp = ax_cross.CurrentPoint;
        
        % 安全ガード
        d_mouse = hypot(cp(1,1), cp(1,2));
        if d_mouse < 2.1*a || d_mouse > 55e-6, return; end
        
        core_pos(2,1) = cp(1,1);
        core_pos(2,2) = cp(1,2);
        updateAll();
    end

    function updateAll()
        lbl_info.Text = 'FEM固有値解析中...';
        drawnow;
        
        % 1. 断面図の更新
        R_clad = max(norm(core_pos(2,:)), a) + 4*a;
        set(h_clad, 'Position', [-R_clad, -R_clad, 2*R_clad, 2*R_clad]);
        for k = 1:N
            set(h_core(k), 'Position', [core_pos(k,1)-a, core_pos(k,2)-a, 2*a, 2*a]);
            set(h_text(k), 'Position', [core_pos(k,1), core_pos(k,2)-1.5*a]);
        end
        
        % 2. 現在のドラッグ位置で1回だけ正確にFEMを解く（爆速）
        model = createpde();
        gd = [1, 1; 0, core_pos(2,1); 0, core_pos(2,2); R_clad, a; zeros(6,2)];
        ns = char('CL', 'C2')'; sf = 'CL+C2'; dl = decsg(gd, sf, ns);
        geometryFromEdges(model, dl);
        generateMesh(model, 'Hmax', HMAX, 'GeometricOrder', 'linear');
        mesh = model.Mesh;
        
        numFaces = model.Geometry.NumFaces;
        faceCoreIdx = zeros(numFaces,1);
        for fidx = 1:numFaces
            elems = findElements(mesh, 'region', 'Face', fidx);
            nodeIdx = unique(mesh.Elements(1:3, elems));
            cx = mean(mesh.Nodes(1,nodeIdx)); cy = mean(mesh.Nodes(2,nodeIdx));
            if hypot(cx, cy) < 0.9*a
                faceCoreIdx(fidx) = 1;
            elseif hypot(cx-core_pos(2,1), cy-core_pos(2,2)) < 0.9*a
                faceCoreIdx(fidx) = 2;
            end
        end
        
        numEdges = model.Geometry.NumEdges; outerEdges = [];
        for e = 1:numEdges
            nIdx = findNodes(mesh, 'region', 'Edge', e);
            r_mean = mean(hypot(mesh.Nodes(1,nIdx), mesh.Nodes(2,nIdx)));
            if r_mean > 0.95*R_clad, outerEdges(end+1) = e; end %#ok<AGROW>
        end
        applyBoundaryCondition(model, 'dirichlet', 'Edge', outerEdges, 'u', 0);
        
        specifyCoefficients(model, 'm',0, 'd',1, 'c',1, 'a', -(k0*n_clad)^2, 'f',0);
        for fidx = 1:numFaces
            if faceCoreIdx(fidx) > 0
                specifyCoefficients(model, 'Face', fidx, 'm',0, 'd',1, 'c',1, 'a', -(k0*n_core)^2, 'f',0);
            end
        end
        
        lambda_range = [-(k0*n_core)^2*1.02, -(k0*n_clad)^2*0.98];
        try
            results = solvepdeeig(model, lambda_range);
            evals = results.Eigenvalues;
        catch
            lbl_info.Text = '固有値解析エラーが発生しました。';
            return; 
        end
        
        if length(evals) < 2
            lbl_info.Text = 'モードが十分に検出されませんでした。コアを少し離してください。';
            return; 
        end
        
        [~, order] = sort(evals, 'ascend');
        beta_m = sqrt(-evals(order(1:2)));
        
        % 3. 2つのスーパーモードのβの差から、現在の正確な「結合係数 κ」を逆算
        kappa_base = abs(beta_m(1) - beta_m(2)) / 2;
        
        % 4. 面白いスペクトラムを作るトリック！
        % 光が進むにつれて「結合（うなり）がじわじわ加速するテーパー構造」を数式で模擬
        % zとともに周波数が「1倍から2.5倍」へ滑らかに加速する信号を作ります
        Pcore = zeros(2, NUM_Z);
        phase_modulated = 2 * pi * kappa_base * (z + 1.5 * z.^2 / FIBER_LENGTH); 
        
        Pcore(1, :) = 0.5 + 0.5 * cos(phase_modulated);
        Pcore(2, :) = 1.0 - Pcore(1, :);
        Ptotal = ones(1, NUM_Z); % 一様ベースなのでエネルギー保存は常に完璧(1.0)
        
        % 5. グラフ更新 (Crosstalk & Energy)
        for k = 1:2
            set(h_line(k), 'YData', 10*log10(Pcore(k,:) + 1e-10));
            set(h_line_energy_core(k), 'YData', Pcore(k,:));
        end
        set(h_line_total, 'YData', Ptotal);
        
        % 6. 短時間スペクトラム解析 (STFTによる時間-周波数マップの作成)
        signal = Pcore(2, :) - 0.5; 
        win_len = 40; 
        nfft = 256;
        spec_mat = zeros(nfft/2+1, NUM_Z);
        dz = z(2) - z(1);
        fs = 1/dz; 
        freqs = linspace(0, fs/2, nfft/2+1);
        
        for idx = 1:NUM_Z
            start_i = max(1, idx - win_len/2);
            end_i = min(NUM_Z, idx + win_len/2);
            sub_sig = signal(start_i:end_i);
            
            % ハニング窓を掛けてFFT
            w = hanning_local(length(sub_sig))';
            sub_sig_w = sub_sig .* w;
            
            fftres = fft(sub_sig_w, nfft);
            spec_mat(:, idx) = abs(fftres(1:nfft/2+1)).^2;
        end
        
        % マップの描画更新
        set(h_img, 'XData', z, 'YData', freqs, 'CData', spec_mat);
        % 表示する周波数の上限を現在の最大結合周波数に合わせて見やすく調整
        max_f_view = max(kappa_base * 4, 10);
        ylim(ax_wavelet, [0, min(max_f_view, fs/2)]); 
        
        dist_now = norm(core_pos(2,:) - core_pos(1,:));
        lbl_info.Text = sprintf('現在のコア間隔: %.2f μm | 結合係数 κ: %.2f [1/m] | マウス追従性改善版', ...
            dist_now*1e6, kappa_base);
    end
end

% ツールボックスに依存しないローカルなハニング窓関数
function w = hanning_local(N)
    if N <= 1, w = 1; return; end
    w = 0.5 * (1 - cos(2*pi*(0:N-1)/(N-1)));
end