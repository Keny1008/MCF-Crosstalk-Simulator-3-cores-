% グラフ描画の完成版イメージ
figure('Position', [100, 100, 900, 600]); 
for core_idx = 2:NUM_CORES
    subplot(2, 3, core_idx - 1); 
    
    % 1. シミュレーション結果（青い線）
    semilogx(bend_diameters_mm, crosstalk_results(:, core_idx), '-b', 'LineWidth', 1.5);
    hold on;
    
    % 2. 論文の実測データ（黒い点）※Core 2のサンプル値
    % 実際には論文のFig.9から値を読み取った配列を入れてください
    measured_x = [100, 200, 300, 500, 700, 1000, 2000]; 
    measured_y = [-48, -45, -43, -41, -54, -57, -59]; 
    semilogx(measured_x, measured_y, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
    
    grid on;
    title(sprintf('Core %d', core_idx));
    xlabel('Bending diameter [mm]'); ylabel('Crosstalk [dB]');
    xlim([100, 2000]); ylim([-60, -20]);
    xticks([100, 200, 500, 1000, 2000]); 
end