function run_mcf_sim()
    % --- 物理パラメータ ---
    global PI WAVELENGTH BETA N_CORE N_CLAD CORE_DIA FIBER_LENGTH
    PI = 3.141592653589793;
    WAVELENGTH = 1550e-9; BETA = 2 * PI / WAVELENGTH;
    N_CLAD = 1.50; DELTA = 0.004; N_CORE = N_CLAD * (1.0 + DELTA);
    CORE_DIA = 8.0e-6; FIBER_LENGTH = 100e-6; 
    
    % 初期状態
    pos1 = -20e-6; pos2 = 20e-6; % コア中心のx座標
    
    % 図の作成
    fig = figure('Name', 'MCF Interactive Simulator', 'NumberTitle', 'off');
    
    % マウス操作の設定
    set(fig, 'WindowButtonDownFcn', @(src,event) start_drag(src));
    set(fig, 'WindowButtonUpFcn', @(src,event) stop_drag(src));
    
    % ドラッグ用変数
    handles.dragging = 0;
    handles.pos2 = pos2;
    guidata(fig, handles);
    
    update_plot(fig, pos1, pos2);
end

function start_drag(fig)
    cp = get(gca, 'CurrentPoint');
    h = guidata(fig);
    if abs(cp(1,1) - h.pos2) < 10e-6, h.dragging = 1; guidata(fig, h); end
end

function stop_drag(fig)
    h = guidata(fig); h.dragging = 0; guidata(fig, h);
end

function drag_core(fig)
    h = guidata(fig);
    if h.dragging
        cp = get(gca, 'CurrentPoint');
        h.pos2 = cp(1,1);
        guidata(fig, h);
        update_plot(fig, -20e-6, h.pos2);
    end
end

% 画面を更新する関数
function update_plot(fig, pos1, pos2)
    global PI WAVELENGTH BETA N_CORE N_CLAD CORE_DIA FIBER_LENGTH
    
    % パラメータ計算
    pitch = abs(pos2 - pos1);
    V = (2*PI/WAVELENGTH)*(CORE_DIA/2)*sqrt(N_CORE^2 - N_CLAD^2);
    u = fzero(@(u_val) u_val*besselj(1,u_val)/besselj(0,u_val) - sqrt(V^2-u_val^2)*besselk(1,sqrt(V^2-u_val^2))/besselk(0,sqrt(V^2-u_val^2)), 2.4);
    v = sqrt(V^2 - u^2); NA = sqrt(N_CORE^2 - N_CLAD^2);
    kappa = (NA/((CORE_DIA/2)*N_CORE)) * (u^2/V^3) * (besselk(0, v*pitch/(CORE_DIA/2)) / besselk(1, v)^2);
    
    z = linspace(0, FIBER_LENGTH, 1000);
    A1 = cos(kappa*z); A2 = sin(kappa*z);
    
    % グラフ描画
    subplot(3,1,1); plot(z*1e6, A1.*sin(BETA*z), 'b', z*1e6, A2.*sin(BETA*z), 'r');
    title(['Electric Field - Pitch: ', num2str(pitch*1e6, '%.1f'), 'um']);
    subplot(3,1,2); plot(z*1e6, A1.^2, 'b', z*1e6, A2.^2, 'r'); title('Power');
    subplot(3,1,3); plot(z*1e6, 10*log10(max(1e-10, A2.^2./(A1.^2+1e-10)))); title('XT [dB]'); ylim([-60 0]);
end