function [bestTbl, rankTbl] = build_strategy_ranking(allSummary, outputDir)
% Build per-group and global ranking from summary table.

if isempty(allSummary)
    bestTbl = table();
    rankTbl = table();
    return;
end

mustCols = {'recipe','lot','parity','signal','strategy','rmse_global','corr_global','bias_rms','pv_mean','rel_pv_change_mean','cv_mean'};
for i = 1:numel(mustCols)
    if ~ismember(mustCols{i}, allSummary.Properties.VariableNames)
        error('Missing required column: %s', mustCols{i});
    end
end

[groupID, recipe, lot, parity, signal] = findgroups(allSummary.recipe, allSummary.lot, allSummary.parity, allSummary.signal);
ug = unique(groupID);
rows = cell(numel(ug), 7);

for i = 1:numel(ug)
    g = ug(i);
    idx = groupID == g;
    sub = allSummary(idx, :);
    [bestRmse, k] = min(sub.rmse_global);
    rows(i,:) = {recipe(i), lot(i), parity(i), signal(i), sub.strategy{k}, bestRmse, sub.corr_global(k)};
end

bestTbl = cell2table(rows, 'VariableNames', ...
    {'recipe','lot','parity','signal','best_strategy','best_rmse','best_corr'});

rankTbl = groupsummary(allSummary, 'strategy', 'mean', {'rmse_global','corr_global','bias_rms','pv_mean','rel_pv_change_mean','cv_mean'});
rankTbl = sortrows(rankTbl, {'mean_rmse_global','mean_cv_mean','mean_bias_rms'}, {'ascend','ascend','ascend'});
rankTbl.rank = (1:height(rankTbl))';
rankTbl = movevars(rankTbl, 'rank', 'Before', 1);

writetable(bestTbl, fullfile(outputDir, 'best_strategy_overall.csv'));
writetable(rankTbl, fullfile(outputDir, 'strategy_rank_overall.csv'));

fid = fopen(fullfile(outputDir, 'strategy_report.txt'), 'w');
fprintf(fid, 'Strategy ranking report\n');
fprintf(fid, '======================\n\n');
if ~isempty(rankTbl)
    fprintf(fid, 'Global best strategy (by mean RMSE): %s\n', rankTbl.strategy{1});
    fprintf(fid, 'Mean RMSE = %.6g, Mean Corr = %.6g\n\n', rankTbl.mean_rmse_global(1), rankTbl.mean_corr_global(1));
end
fprintf(fid, 'Detailed ranking:\n');
for i = 1:height(rankTbl)
    fprintf(fid, '%d) %s | RMSE=%.6g Corr=%.6g BiasRMS=%.6g CV=%.6g\n', ...
        rankTbl.rank(i), rankTbl.strategy{i}, rankTbl.mean_rmse_global(i), rankTbl.mean_corr_global(i), ...
        rankTbl.mean_bias_rms(i), rankTbl.mean_cv_mean(i));
end
fclose(fid);
end
