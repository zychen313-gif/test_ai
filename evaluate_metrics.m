function M = evaluate_metrics(errCell, predCell, actCell, rCell, tCell, cfg)
err = vertcat(errCell{:});
pred = vertcat(predCell{:});
act = vertcat(actCell{:});

M.rmse_global = sqrt(mean(err.^2));
if std(pred) > 0 && std(act) > 0
    C = corrcoef(pred, act);
    M.corr_global = C(1,2);
else
    M.corr_global = NaN;
end
M.bias_rms = sqrt(mean((mean(pred)-mean(act)).^2));

pvAct = cellfun(@(x) max(x)-min(x), actCell);
pvErr = cellfun(@(x) max(x)-min(x), errCell);
M.pv_mean = mean(pvErr);
M.rel_pv_change_mean = mean(pvErr ./ max(pvAct, eps));

[rMap, tMap, stdMap, cvMap, radialStd] = build_spatial_maps(errCell, actCell, rCell, tCell, cfg);
M.rMap = rMap; M.tMap = tMap;
M.stdMap = stdMap;
M.cvMap = cvMap;
M.radialStd = radialStd;
M.cv_mean = mean(cvMap(isfinite(cvMap)), 'omitnan');

M.hist_edges = linspace(min(err), max(err), 40);
M.hist_counts = histcounts(err, M.hist_edges);
M.err_all = err;
end

function [rGrid, tGrid, stdMap, cvMap, radialStd] = build_spatial_maps(errCell, actCell, rCell, tCell, cfg)
rGrid = linspace(cfg.rMin, cfg.rMax, 60);
tGrid = linspace(0, 360, 120);
stdMap = nan(numel(rGrid), numel(tGrid));
cvMap = nan(numel(rGrid), numel(tGrid));

for ir = 1:numel(rGrid)-1
    for it = 1:numel(tGrid)-1
        ebin = [];
        abin = [];
        for k = 1:numel(errCell)
            in = rCell{k} >= rGrid(ir) & rCell{k} < rGrid(ir+1) & tCell{k} >= tGrid(it) & tCell{k} < tGrid(it+1);
            ebin = [ebin; errCell{k}(in)]; %#ok<AGROW>
            abin = [abin; actCell{k}(in)]; %#ok<AGROW>
        end
        if numel(ebin) >= 5
            stdMap(ir,it) = std(ebin, 'omitnan');
            cvMap(ir,it) = std(ebin, 'omitnan') / max(abs(mean(abin, 'omitnan')), eps);
        end
    end
end
radialStd = mean(stdMap,2,'omitnan');
end
