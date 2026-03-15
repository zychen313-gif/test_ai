function save_strategy_outputs(M, strategy, resultDir)
outCsv = fullfile(resultDir, sprintf('%s_metrics.csv', strategy));
T = table(M.rmse_global, M.corr_global, M.bias_rms, M.pv_mean, M.rel_pv_change_mean, M.cv_mean, ...
    'VariableNames', {'rmse_global','corr_global','bias_rms','pv_mean','rel_pv_change_mean','cv_mean'});
writetable(T, outCsv);

fig1 = figure('Visible','off');
imagesc(M.tMap, M.rMap, M.stdMap); axis xy; colorbar;
xlabel('\theta (deg)'); ylabel('r (mm)'); title([strategy ' std-map']);
saveas(fig1, fullfile(resultDir, sprintf('%s_std_map.png', strategy))); close(fig1);

fig2 = figure('Visible','off');
imagesc(M.tMap, M.rMap, M.cvMap); axis xy; colorbar;
xlabel('\theta (deg)'); ylabel('r (mm)'); title([strategy ' CV-map']);
saveas(fig2, fullfile(resultDir, sprintf('%s_cv_map.png', strategy))); close(fig2);

fig3 = figure('Visible','off');
bar(M.hist_edges(1:end-1), M.hist_counts, 'histc');
xlabel('Prediction error'); ylabel('Count'); title([strategy ' error histogram']);
saveas(fig3, fullfile(resultDir, sprintf('%s_err_hist.png', strategy))); close(fig3);

fig4 = figure('Visible','off');
plot(M.rMap, M.radialStd, 'LineWidth', 1.5);
xlabel('r (mm)'); ylabel('std(error)'); title([strategy ' radial consistency profile']); grid on;
saveas(fig4, fullfile(resultDir, sprintf('%s_radial_profile.png', strategy))); close(fig4);
end
