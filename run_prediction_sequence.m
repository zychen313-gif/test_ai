function [summaryTbl, metricsByStrategy] = run_prediction_sequence(dataSeq, waferNames, signalType, cfg, resultDir)
strategies = cfg.strategies;
metricsByStrategy = struct();
summaryRows = {};

% Batch average shape (edge-only) for current parity/signal
save_average_shape(dataSeq, signalType, cfg, resultDir);


for s = 1:numel(strategies)
    strategy = strategies{s};
    predErrAll = {};
    predValAll = {};
    actValAll = {};
    rAll = {};
    thAll = {};
    waferUsed = {};

    for i = 2:numel(dataSeq)
        train = dataSeq(1:i-1);
        test = dataSeq{i};

        [xE, yE, rE, thE, zTest] = get_edge_signal(test, signalType);
        if numel(zTest) < 50
            continue;
        end
        model = fit_edge_model(train, signalType, strategy, cfg);
        zPred = predict_edge_model(model, xE, yE, rE, thE, strategy);

        err = zPred - zTest;
        predErrAll{end+1} = err; %#ok<AGROW>
        predValAll{end+1} = zPred; %#ok<AGROW>
        actValAll{end+1} = zTest; %#ok<AGROW>
        rAll{end+1} = rE; %#ok<AGROW>
        thAll{end+1} = thE; %#ok<AGROW>
        waferUsed{end+1} = waferNames{i}; %#ok<AGROW>
    end

    if isempty(predErrAll)
        continue;
    end

    M = evaluate_metrics(predErrAll, predValAll, actValAll, rAll, thAll, cfg);
    metricsByStrategy.(strategy) = M;
    save_strategy_outputs(M, strategy, resultDir);

    summaryRows(end+1,:) = {strategy, M.rmse_global, M.corr_global, M.bias_rms, M.pv_mean, M.rel_pv_change_mean, M.cv_mean}; %#ok<AGROW>
end

if isempty(summaryRows)
    summaryTbl = table();
else
    summaryTbl = cell2table(summaryRows, 'VariableNames', ...
        {'strategy','rmse_global','corr_global','bias_rms','pv_mean','rel_pv_change_mean','cv_mean'});
end

if ~isempty(summaryTbl)
    writetable(summaryTbl, fullfile(resultDir, 'summary.csv'));
end
end

function [xE, yE, rE, thE, zE] = get_edge_signal(d, signalType)
mask = d.edgeMask;
xE = d.x(mask); yE = d.y(mask); rE = d.r(mask); thE = d.theta(mask);
switch signalType
    case 'X', zE = d.zX(mask);
    case 'Y', zE = d.zY(mask);
    otherwise, error('Unknown signal type');
end
end


function save_average_shape(dataSeq, signalType, cfg, resultDir)
rGrid = linspace(cfg.rMin, cfg.rMax, 80);
tGrid = linspace(0, 360, 180);
acc = zeros(numel(rGrid)-1, numel(tGrid)-1);
cnt = zeros(numel(rGrid)-1, numel(tGrid)-1);
for i = 1:numel(dataSeq)
    d = dataSeq{i};
    m = d.edgeMask;
    r = d.r(m); t = d.theta(m);
    if signalType == 'X'
        z = d.zX(m);
    else
        z = d.zY(m);
    end
    for ir = 1:numel(rGrid)-1
        for it = 1:numel(tGrid)-1
            in = r>=rGrid(ir) & r<rGrid(ir+1) & t>=tGrid(it) & t<tGrid(it+1);
            if any(in)
                acc(ir,it) = acc(ir,it) + mean(z(in),'omitnan');
                cnt(ir,it) = cnt(ir,it) + 1;
            end
        end
    end
end
meanMap = acc ./ max(cnt,1);
fig = figure('Visible','off');
imagesc(tGrid(1:end-1), rGrid(1:end-1), meanMap); axis xy; colorbar;
xlabel('\theta (deg)'); ylabel('r (mm)'); title(['Average edge shape (' signalType ')']);
saveas(fig, fullfile(resultDir, sprintf('average_shape_%s.png', signalType))); close(fig);

writematrix(meanMap, fullfile(resultDir, sprintf('average_shape_%s.csv', signalType)));
end
