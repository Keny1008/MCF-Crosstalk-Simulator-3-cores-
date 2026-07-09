function drag_interactive_mcf_crosstalk_FEM()
    % --- 1. 固定の物理パラメータ設定 ---
    e0 = 8.85418781e-12;
    m0 = 1.25663706e-6;
    c_speed = 1.0 / sqrt(e0 * m0);
    f = 193e12;
    lambda = c_speed / f;
    k0 = 2*pi/lambda;

    n_clad = input('クラッドの屈折率を入力: ');
    n_core = input('コアの屈折率を入力: ');
    a = input('コアの半径 [μm]を入力: ') * 1e-6;
    N = input('コア数を入力: ');

    FIBER_LENGTH = 10000;
    NUM_Z = 400;                 % FEM版は重いのでzは粗めに(必要なら増やす)
    z = linspace(0, FIBER_LENGTH, NUM_Z);

    HMAX = a / 2.2;              % メッシュの粗さ。小さくすると精度↑・速度↓

    % --- 2. GUIの作成 ---
    fig = uifigure('Name', 'FEM-based Draggable MCF Crosstalk (N-core)', 'Position', [100, 100, 1150, 550]);
    gl = uigridlayout(fig, [2, 2]);
    gl.RowHeight = {'1x', 40};
    gl.ColumnWidth = {'1.5x', '1x'};

    ax_graph = uiaxes(gl);
    ax_graph.Layout.Row = 1; ax_graph.Layout.Column = 1;
    xlabel(ax_graph, 'Distance [m]'); ylabel(ax_graph, 'Crosstalk [dB]');
    grid(ax_graph, 'on'); hold(ax_graph, 'on');

    ax_cross = uiaxes(gl);
    ax_cross.Layout.Row = 1; ax_cross.Layout.Column = 2;
    xlabel(ax_cross, 'x [m]'); ylabel(ax_cross, 'y [m]');
    axis(ax_cross, 'equal'); grid(ax_cross, 'on'); hold(ax_cross, 'on');
    title(ax_cross, 'コアをドラッグ (FEMで断面全体を再計算)');
    xlim(ax_cross, [-100e-6, 140e-6]); ylim(ax_cross, [-100e-6, 160e-6]);

    lbl_info = uilabel(gl, 'Text', '初期化中...', 'FontWeight', 'bold', 'FontSize', 11, 'WordWrap', 'on');
    lbl_info.Layout.Row = 2; lbl_info.Layout.Column = [1 2];

    % --- 3. コア座標の初期設定(リング状配置、N=3なら正三角形と同等) ---
    p_init = 39.2e-6;
    core_pos = zeros(N, 2);
    core_pos(1,:) = [0, 0];
    if N > 1
        theta = linspace(0, 2*pi, N)'; theta(end) = [];
        core_pos(2:end,1) = p_init * cos(theta);
        core_pos(2:end,2) = p_init * sin(theta);
    end

    % --- 4. 描画オブジェクト(断面図・グラフ)の初期化 ---
    h_clad = rectangle(ax_cross, 'Position', [0 0 1 1], 'Curvature', [1 1], ...
                        'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'k', 'PickableParts', 'none');
    colors = lines(N);
    h_core = gobjects(N,1); h_text = gobjects(N,1);
    for k = 1:N
        col = colors(k,:); if k==1, col = [1 0 0]; end
        h_core(k) = rectangle(ax_cross, 'Position', [0 0 1 1], 'Curvature', [1 1], ...
                               'FaceColor', col, 'EdgeColor', 'k', 'PickableParts', 'all');
        h_text(k) = text(ax_cross, 0, 0, sprintf('Core %d', k), ...
                          'HorizontalAlignment', 'center', 'PickableParts', 'none', 'FontWeight', 'bold');
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

    % --- 5. ドラッグ操作のコールバック ---
    current_drag_core = 0;
    for k = 1:N
        h_core(k).ButtonDownFcn = @(~,~) startDrag(k);
    end
    fig.WindowButtonMotionFcn = @dragging;
    fig.WindowButtonUpFcn = @stopDrag;

    updateAll();

    function startDrag(core_idx), current_drag_core = core_idx; end
    function stopDrag(~,~), current_drag_core = 0; end
    function dragging(~,~)
        if current_drag_core == 0, return; end
        cp = ax_cross.CurrentPoint;
        core_pos(current_drag_core,1) = cp(1,1);
        core_pos(current_drag_core,2) = cp(1,2);
        updateAll();
    end

    % --- 6. メインの更新関数(ジオメトリ生成→FEM固有値解析→伝搬計算) ---
    function updateAll()
        lbl_info.Text = 'FEM計算中...';
        drawnow;

        model = createpde();   % geometryFromEdgesは同一modelに2回呼べないため毎回新規作成

        % --- 6-1. 断面図(円の位置)の更新 ---
        cent = mean(core_pos,1);
        maxd = max(vecnorm(core_pos - cent, 2, 2));
        R_clad = maxd + 4*a;
        set(h_clad, 'Position', [cent(1)-R_clad, cent(2)-R_clad, 2*R_clad, 2*R_clad]);
        for k = 1:N
            set(h_core(k), 'Position', [core_pos(k,1)-a, core_pos(k,2)-a, 2*a, 2*a]);
            set(h_text(k), 'Position', [core_pos(k,1), core_pos(k,2)-1.5*a]);
        end

        % --- 6-2. decsgでジオメトリ作成(クラッド円 + コア円 x N) ---
        names = [{'CL'}, arrayfun(@(k) sprintf('C%d',k), 1:N, 'UniformOutput', false)];
        gd = zeros(10, N+1);
        gd(:,1) = [1; cent(1); cent(2); R_clad; zeros(6,1)];
        for k = 1:N
            gd(:,k+1) = [1; core_pos(k,1); core_pos(k,2); a; zeros(6,1)];
        end
        ns = char(names)';
        sf = strjoin(names, '+');
        dl = decsg(gd, sf, ns);

        geometryFromEdges(model, dl);
        generateMesh(model, 'Hmax', HMAX, 'GeometricOrder', 'linear');
        mesh = model.Mesh;

        % --- 6-3. 各Faceがどのコアに属するか、中心座標で判定 ---
        numFaces = model.Geometry.NumFaces;
        faceCoreIdx = zeros(numFaces,1);   % 0 = クラッド, k = コアk
        for fidx = 1:numFaces
            elems = findElements(mesh, 'region', 'Face', fidx);
            nodeIdx = unique(mesh.Elements(1:3, elems));
            cx = mean(mesh.Nodes(1,nodeIdx)); cy = mean(mesh.Nodes(2,nodeIdx));
            for k = 1:N
                if hypot(cx-core_pos(k,1), cy-core_pos(k,2)) < 0.9*a
                    faceCoreIdx(fidx) = k; break;
                end
            end
        end

        % --- 6-4. 外周(ディリクレ境界)のエッジを、半径が最大に近いもので判定 ---
        numEdges = model.Geometry.NumEdges;
        outerEdges = [];
        for e = 1:numEdges
            nIdx = findNodes(mesh, 'region', 'Edge', e);
            r_mean = mean(hypot(mesh.Nodes(1,nIdx)-cent(1), mesh.Nodes(2,nIdx)-cent(2)));
            if r_mean > 0.97*R_clad
                outerEdges(end+1) = e; %#ok<AGROW>
            end
        end
        applyBoundaryCondition(model, 'dirichlet', 'Edge', outerEdges, 'u', 0);

        % --- 6-5. 固有値問題の係数設定: -Δψ - k0^2 n(x,y)^2 ψ = λψ  (λ=-β^2) ---
        specifyCoefficients(model, 'm',0, 'd',1, 'c',1, 'a', -(k0*n_clad)^2, 'f',0); % 既定=クラッド
        for fidx = 1:numFaces
            if faceCoreIdx(fidx) > 0
                specifyCoefficients(model, 'Face', fidx, 'm',0, 'd',1, 'c',1, 'a', -(k0*n_core)^2, 'f',0);
            end
        end

        % --- 6-6. 固有値解析(ガイドモード帯 k0*n_clad < beta < k0*n_core 相当) ---
        lambda_range = [-(k0*n_core)^2*1.02, -(k0*n_clad)^2*0.98];
        try
            results = solvepdeeig(model, lambda_range);
        catch ME
            lbl_info.Text = ['固有値解析エラー: ', ME.message];
            return;
        end
        evals = results.Eigenvalues;
        evecs = results.Eigenvectors;

        if length(evals) < N
            lbl_info.Text = sprintf('警告: 見つかったモード数(%d)がコア数(%d)より少ないです。コア間隔やHMAXを調整してください。', length(evals), N);
        end
        Nuse = min(N, length(evals));
        [~, order] = sort(evals, 'ascend');   % 最も負(=beta最大)から
        idx_use = order(1:Nuse);
        beta_m = sqrt(-evals(idx_use));
        Psi = evecs(:, idx_use);

        % --- 6-7. 全域質量行列M(規格化・射影用) ---
        specifyCoefficients(model, 'm',0, 'd',1, 'c',1, 'a',0, 'f',0);
        FEMfull = assembleFEMatrices(model, 'M');
        Mmat = FEMfull.M;

        normM = sqrt(diag(Psi.' * Mmat * Psi));
        Psi = Psi ./ normM.';

        % --- 6-8. 入力励振(Core1中心のガウシアン近似)を射影 ---
        nodesXY = mesh.Nodes;
        w0 = a;
        psi_in = exp(-((nodesXY(1,:)-core_pos(1,1)).^2 + (nodesXY(2,:)-core_pos(1,2)).^2) / w0^2).';
        c_m = Psi.' * Mmat * psi_in;

        % --- 6-9. 各コア領域ごとの質量行列でパワーを算出しながらz方向に伝搬 ---
        phase = exp(-1j * beta_m * z);      % (Nuse x NUM_Z)
        Uz = Psi * (c_m .* phase);          % (Nnodes x NUM_Z)

        % 注意: PDE Toolboxは「同じ係数があるFaceでは0、別のFaceでは非0」という
        % 設定を許可しない(Inconsistent coefficients across subdomains エラー)。
        % そのため、対象外の面はd=0ではなくd=EPS(無視できるほど小さい値)にする。
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
            Mk = FEMk.M;
            Pcore(k,:) = real(sum(conj(Uz) .* (Mk * Uz), 1));
        end
        Pcore = Pcore / max(Pcore(1,1), eps);   % z=0でCore1=1に規格化

        for k = 1:N
            set(h_line(k), 'YData', 10*log10(Pcore(k,:) + 1e-20));
        end

        % --- 6-10. ピッチ情報の表示 ---
        info_str = '';
        for ii = 1:N
            for jj = ii+1:N
                d_ij = norm(core_pos(ii,:) - core_pos(jj,:));
                info_str = [info_str, sprintf('[C%d-C%d]: %.2f \\mum   ', ii, jj, d_ij*1e6)]; %#ok<AGROW>
            end
        end
        lbl_info.Text = ['ピッチ  ', info_str];
    end
end