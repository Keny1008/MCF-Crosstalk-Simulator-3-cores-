function MCF_crosstalk_FEM_energy_check()

    %% --- 1. 物理定数・パラメータ入力(前回コードと同じ体系) ---
    e0 = 8.85418781e-12;
    m0 = 1.25663706e-6;
    c_speed = 1.0 / sqrt(e0 * m0);
    f = 193e12;
    lambda = c_speed / f;
    k0 = 2*pi/lambda;

    n_clad = input('クラッドの屈折率を入力: ');
    n_core = input('コアの屈折率を入力: ');
    N      = input('コア数を入力: ');

    a = zeros(N,1);
    for k = 1:N
        a(k) = input(sprintf('コア%dの半径 [μm]を入力: ', k)) * 1e-6;
    end

    FIBER_LENGTH = 10000;          % [m]
    NUM_Z = 400;
    z = linspace(0, FIBER_LENGTH, NUM_Z);

    R_CLAD = 125e-6 / 2;            % クラッドは直径125μmで固定(半径62.5μm)。以後変更しない。
    HMAX = min(a) / 2.2;             % 最小コアに合わせてメッシュ細かさを決定

    %% --- 2. GUI作成 ---
    fig = uifigure('Name', 'MCF Crosstalk FEM + Energy Conservation Check', ...
                   'Position', [80, 60, 1250, 800]);
    gl = uigridlayout(fig, [3, 2]);
    gl.RowHeight = {'1x', '1x', 40};
    gl.ColumnWidth = {'1.5x', '1x'};

    ax_graph = uiaxes(gl);
    ax_graph.Layout.Row = 1; ax_graph.Layout.Column = 1;
    xlabel(ax_graph, 'Distance [m]'); ylabel(ax_graph, 'Crosstalk [dB]');
    title(ax_graph, '各コアのパワー(クロストーク)');
    grid(ax_graph, 'on'); hold(ax_graph, 'on');

    ax_cross = uiaxes(gl);
    ax_cross.Layout.Row = 1; ax_cross.Layout.Column = 2;
    xlabel(ax_cross, 'x [m]'); ylabel(ax_cross, 'y [m]');
    axis(ax_cross, 'equal'); grid(ax_cross, 'on'); hold(ax_cross, 'on');
    title(ax_cross, 'コアをドラッグ (クラッドは125\mum固定・移動不可)');
    xlim(ax_cross, [-R_CLAD*1.2, R_CLAD*1.2]);
    ylim(ax_cross, [-R_CLAD*1.2, R_CLAD*1.2]);

    ax_energy = uiaxes(gl);
    ax_energy.Layout.Row = 2; ax_energy.Layout.Column = [1 2];
    xlabel(ax_energy, 'Distance [m]'); ylabel(ax_energy, 'Power (normalized)');
    title(ax_energy, 'エネルギー保存チェック: 各コア + 残りクラッド + 合計(=1になるか)');
    grid(ax_energy, 'on'); hold(ax_energy, 'on');
    ylim(ax_energy, [-0.1, 1.2]);

    lbl_info = uilabel(gl, 'Text', '初期化中...', 'FontWeight', 'bold', ...
                        'FontSize', 11, 'WordWrap', 'on');
    lbl_info.Layout.Row = 3; lbl_info.Layout.Column = [1 2];

    %% --- 3. コア初期配置(クラッド内に収まるリング状配置) ---
    p_init = min(0.55 * R_CLAD, 39.2e-6);
    core_pos = zeros(N, 2);
    core_pos(1,:) = [0, 0];
    if N > 1
        theta = linspace(0, 2*pi, N)'; theta(end) = [];
        core_pos(2:end,1) = p_init * cos(theta);
        core_pos(2:end,2) = p_init * sin(theta);
    end
    core_pos_prev = core_pos;   % 不正なドラッグ・計算失敗時のロールバック用

    %% --- 4. 描画オブジェクトの初期化 ---
    % クラッドはドラッグ不可(PickableParts=none) & 半径固定で以後一切更新しない
    rectangle(ax_cross, 'Position', [-R_CLAD, -R_CLAD, 2*R_CLAD, 2*R_CLAD], ...
        'Curvature', [1 1], 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'k', ...
        'PickableParts', 'none');

    colors = lines(N);
    h_core = gobjects(N,1); h_text = gobjects(N,1);
    for k = 1:N
        col = colors(k,:); if k==1, col = [1 0 0]; end
        h_core(k) = rectangle(ax_cross, 'Position', ...
            [core_pos(k,1)-a(k), core_pos(k,2)-a(k), 2*a(k), 2*a(k)], ...
            'Curvature', [1 1], 'FaceColor', col, 'EdgeColor', 'k', 'PickableParts', 'all');
        h_text(k) = text(ax_cross, core_pos(k,1), core_pos(k,2)-1.5*a(k), ...
            sprintf('Core %d', k), 'HorizontalAlignment', 'center', ...
            'PickableParts', 'none', 'FontWeight', 'bold');
    end

    h_line = gobjects(N,1); leg_labels = cell(N,1);
    for k = 1:N
        col = colors(k,:); lw = 1.5;
        if k == 1, col = [0 0 0]; lw = 2; end
        h_line(k) = plot(ax_graph, z, nan(size(z)), 'Color', col, 'LineWidth', lw);
        if k == 1, leg_labels{k} = 'Core 1 (Input)'; else, leg_labels{k} = sprintf('Core %d', k); end
    end
    ylim(ax_graph, [-160, 5]); xlim(ax_graph, [0, FIBER_LENGTH]);
    legend(ax_graph, h_line, leg_labels, 'Location', 'southwest');

    h_eline = gobjects(N,1);
    for k = 1:N
        col = colors(k,:); if k==1, col = [1 0 0]; end
        h_eline(k) = plot(ax_energy, z, nan(size(z)), 'Color', col, 'LineWidth', 1.3);
    end
    h_eclad  = plot(ax_energy, z, nan(size(z)), 'Color', [0.4 0.4 0.4], 'LineWidth', 1.3, 'LineStyle', '--');
    h_etotal = plot(ax_energy, z, nan(size(z)), 'Color', [0 0.6 0], 'LineWidth', 2.2, 'LineStyle', ':');
    plot(ax_energy, z, ones(size(z)), 'Color', [1 0 1], 'LineWidth', 0.8, 'LineStyle', '-.'); % 基準線(=1)

    leg2 = [h_eline; h_eclad; h_etotal];
    leg2_labels = [arrayfun(@(k) sprintf('Core %d power', k), 1:N, 'UniformOutput', false), ...
                   {'Cladding remainder'}, {'Sum (core+clad, should = 1)'}];
    legend(ax_energy, leg2, leg2_labels, 'Location', 'eastoutside');

    %% --- 5. ドラッグ操作 ---
    current_drag_core = 0;
    for k = 1:N
        h_core(k).ButtonDownFcn = @(~,~) startDrag(k);
    end
    fig.WindowButtonMotionFcn = @dragging;
    fig.WindowButtonUpFcn = @stopDrag;

    updateAll();

    function startDrag(core_idx)
        current_drag_core = core_idx;
        core_pos_prev = core_pos;
    end
    function stopDrag(~,~)
        current_drag_core = 0;
    end
    function dragging(~,~)
        if current_drag_core == 0, return; end
        cp = ax_cross.CurrentPoint;
        k = current_drag_core;
        newx = cp(1,1); newy = cp(1,2);

        % --- クラッドは固定のまま、コアがクラッド内に収まるようクランプ ---
        r_new = hypot(newx, newy);
        r_limit = R_CLAD - a(k) - 1e-9;
        if r_limit <= 0
            return;   % コア半径がクラッドより大きい異常設定
        end
        if r_new > r_limit
            scale = r_limit / max(r_new, eps);
            newx = newx * scale; newy = newy * scale;
        end

        % --- 他コアとの重なりを防止(重なる移動は無視) ---
        for j = 1:N
            if j == k, continue; end
            if hypot(newx-core_pos(j,1), newy-core_pos(j,2)) < (a(k)+a(j))*1.05
                return;
            end
        end

        core_pos(k,:) = [newx, newy];
        updateAll();
    end

    %% --- 6. メイン計算(ジオメトリ生成 → FEM固有値解析 → 伝搬 → エネルギー検証) ---
    function updateAll()
        lbl_info.Text = 'FEM計算中...';
        drawnow;
        model = createpde();

        % --- 6-1. 断面図の更新(クラッドは触らない、コアのみ更新) ---
        for k = 1:N
            set(h_core(k), 'Position', [core_pos(k,1)-a(k), core_pos(k,2)-a(k), 2*a(k), 2*a(k)]);
            set(h_text(k), 'Position', [core_pos(k,1), core_pos(k,2)-1.5*a(k)]);
        end

        % --- 6-2. decsgでジオメトリ作成(クラッド円[固定,半径125/2um] + コア円xN) ---
        names = [{'CL'}, arrayfun(@(k) sprintf('C%d',k), 1:N, 'UniformOutput', false)];
        gd = zeros(10, N+1);
        gd(:,1) = [1; 0; 0; R_CLAD; zeros(6,1)];
        for k = 1:N
            gd(:,k+1) = [1; core_pos(k,1); core_pos(k,2); a(k); zeros(6,1)];
        end
        ns = char(names)';
        sf = strjoin(names, '+');
        dl = decsg(gd, sf, ns);
        geometryFromEdges(model, dl);
        generateMesh(model, 'Hmax', HMAX, 'GeometricOrder', 'linear');
        mesh = model.Mesh;

        % --- 6-3. 各Faceがどのコアに属するか判定(0=クラッド) ---
        numFaces = model.Geometry.NumFaces;
        faceCoreIdx = zeros(numFaces,1);
        for fidx = 1:numFaces
            elems = findElements(mesh, 'region', 'Face', fidx);
            nodeIdx = unique(mesh.Elements(1:3, elems));
            cx = mean(mesh.Nodes(1,nodeIdx)); cy = mean(mesh.Nodes(2,nodeIdx));
            for k = 1:N
                if hypot(cx-core_pos(k,1), cy-core_pos(k,2)) < 0.9*a(k)
                    faceCoreIdx(fidx) = k; break;
                end
            end
        end

        % --- 6-4. 外周(ディリクレ境界)判定 ---
        numEdges = model.Geometry.NumEdges;
        outerEdges = [];
        for e = 1:numEdges
            nIdx = findNodes(mesh, 'region', 'Edge', e);
            r_mean = mean(hypot(mesh.Nodes(1,nIdx), mesh.Nodes(2,nIdx)));
            if r_mean > 0.97*R_CLAD
                outerEdges(end+1) = e; %#ok<AGROW>
            end
        end
        applyBoundaryCondition(model, 'dirichlet', 'Edge', outerEdges, 'u', 0);

        % --- 6-5. 固有値問題係数: -Δψ - k0^2 n(x,y)^2 ψ = λψ  (λ = -β^2) ---
        specifyCoefficients(model, 'm',0, 'd',1, 'c',1, 'a', -(k0*n_clad)^2, 'f',0);
        for fidx = 1:numFaces
            if faceCoreIdx(fidx) > 0
                specifyCoefficients(model, 'Face', fidx, 'm',0, 'd',1, 'c',1, ...
                    'a', -(k0*n_core)^2, 'f',0);
            end
        end

        lambda_range = [-(k0*n_core)^2*1.02, -(k0*n_clad)^2*0.98];
        try
            results = solvepdeeig(model, lambda_range);
        catch ME
            lbl_info.Text = ['固有値解析エラー: ', ME.message];
            core_pos = core_pos_prev;
            return;
        end
        evals = results.Eigenvalues;
        evecs = results.Eigenvectors;

        if length(evals) < N
            lbl_info.Text = sprintf(['警告: 見つかったモード数(%d)がコア数(%d)より少ないです。' ...
                'コア間隔・半径・HMAXを調整してください。'], length(evals), N);
        end
        Nuse = min(N, length(evals));
        [~, order] = sort(evals, 'ascend');     % 最も負(=betaが最大)のものから
        idx_use = order(1:Nuse);
        beta_m = sqrt(-evals(idx_use));
        Psi = evecs(:, idx_use);

        % --- 6-6. 全域質量行列(モード規格化・射影・全体パワー検証用) ---
        specifyCoefficients(model, 'm',0, 'd',1, 'c',1, 'a',0, 'f',0);
        FEMfull = assembleFEMatrices(model, 'M');
        Mmat = FEMfull.M;
        normM = sqrt(diag(Psi.' * Mmat * Psi));
        Psi = Psi ./ normM.';   % 各スーパーモードを ∫|psi|^2 dA = 1 に規格化(正規直交化)

        % --- 6-7. 入力励振(Core1中心のガウシアン、入力パワーを1に規格化) ---
        nodesXY = mesh.Nodes;
        w0 = a(1);
        psi_in = exp(-((nodesXY(1,:)-core_pos(1,1)).^2 + (nodesXY(2,:)-core_pos(1,2)).^2) / w0^2).';
        Pin_raw = psi_in.' * Mmat * psi_in;
        psi_in = psi_in / sqrt(Pin_raw);            % ∫|psi_in|^2 dA = 1 に規格化

        c_m = Psi.' * Mmat * psi_in;                % 各スーパーモードへの射影係数
        c_m = c_m / sqrt(sum(abs(c_m).^2));         % 計算対象のNuseモードのみで入力パワー100%と仮定
                                                     % (これにより理論上の合計パワーが厳密に1になる)

        % --- 6-8. z方向伝搬(位相のみ変化する無損失コヒーレント重ね合わせ) ---
        phase = exp(-1j * beta_m * z);
        Uz = Psi * (c_m .* phase);                  % (Nnodes x NUM_Z)

        % --- 6-9. 各コア領域のパワー(質量行列で領域積分) ---
        EPS_D = 1e-12;
        Pcore = zeros(N, NUM_Z);
        for k = 1:N
            specifyCoefficients(model, 'm',0, 'd',EPS_D, 'c',1, 'a',0, 'f',0);
            for fidx = 1:numFaces
                if faceCoreIdx(fidx) == k
                    specifyCoefficients(model, 'Face', fidx, 'm',0, 'd',1, 'c',1, 'a',0, 'f',0);
                end
            end
            FEMk = assembleFEMatrices(model, 'M');
            Pcore(k,:) = real(sum(conj(Uz) .* (FEMk.M * Uz), 1));
        end

        % --- 6-10. 残りクラッド領域のパワー(コア領域とは別に独立にアセンブル) ---
        specifyCoefficients(model, 'm',0, 'd',EPS_D, 'c',1, 'a',0, 'f',0);
        for fidx = 1:numFaces
            if faceCoreIdx(fidx) == 0
                specifyCoefficients(model, 'Face', fidx, 'm',0, 'd',1, 'c',1, 'a',0, 'f',0);
            end
        end
        FEMclad = assembleFEMatrices(model, 'M');
        Pclad = real(sum(conj(Uz) .* (FEMclad.M * Uz), 1));

        % --- 6-11. 全域パワー(全域Mで直接計算=理論上つねに1)との比較 ---
        Ptotal_direct = real(sum(conj(Uz) .* (Mmat * Uz), 1));   % 理論上 z によらず 1
        Ptotal_sum    = sum(Pcore, 1) + Pclad;                    % 各領域を個別に積分した和(独立計算)

        % --- 6-12. クロストークグラフ(dB, Core1のz=0でのパワーを基準) ---
        Pcore_norm = Pcore ./ max(Pcore(1,1), eps);
        for k = 1:N
            set(h_line(k), 'YData', 10*log10(Pcore_norm(k,:) + 1e-20));
        end

        % --- 6-13. エネルギー保存グラフ ---
        for k = 1:N
            set(h_eline(k), 'YData', Pcore(k,:));
        end
        set(h_eclad,  'YData', Pclad);
        set(h_etotal, 'YData', Ptotal_sum);

        max_diff = max(abs(Ptotal_sum - Ptotal_direct));   % 独立計算どうしのズレ=実装の自己無矛盾性の指標

        % --- 6-14. ピッチ・検証情報の表示 ---
        info_str = '';
        for ii = 1:N
            for jj = ii+1:N
                d_ij = norm(core_pos(ii,:) - core_pos(jj,:));
                info_str = [info_str, sprintf('[C%d-C%d]: %.2f \\mum   ', ii, jj, d_ij*1e6)]; %#ok<AGROW>
            end
        end
        lbl_info.Text = sprintf(['ピッチ  %s | エネルギー保存チェック: 独立に足し合わせた合計パワー ' ...
            'と全域直接計算値との最大誤差 = %.3e (理想値0。数値誤差程度に収まっていれば実装は自己無矛盾)'], ...
            info_str, max_diff);
    end
end