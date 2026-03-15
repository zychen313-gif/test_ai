function summaryTbl = process_recipe(recipePath, outputDir, cfg)
recipeName = get_last_path_name(recipePath);
lotDirs = dir(recipePath);
lotDirs = lotDirs([lotDirs.isdir]);
lotDirs = lotDirs(~ismember({lotDirs.name},{'.','..'}));
summaryTbl = table();

for l = 1:numel(lotDirs)
    lotPath = fullfile(recipePath, lotDirs(l).name);
    lotName = lotDirs(l).name;
    wafers = collect_wafer_files(lotPath);
    if numel(wafers) < 2
        continue;
    end

    data = cell(numel(wafers), 1);
    for k = 1:numel(wafers)
        data{k} = load_and_preprocess_wafer(wafers{k}, cfg);
    end

    idxOdd = 1:2:numel(wafers);
    idxEven = 2:2:numel(wafers);

    for parityCell = {'odd','even'}
        parity = parityCell{1};
        switch parity
            case 'odd', idx = idxOdd;
            case 'even', idx = idxEven;
        end
        if numel(idx) < 2
            continue;
        end
        subset = data(idx);
        waferNames = wafers(idx);
        for signalCell = {'X','Y'}
            signalType = signalCell{1};
            resultDir = fullfile(outputDir, recipeName, lotName, parity, signalType);
            if ~exist(resultDir, 'dir')
                mkdir(resultDir);
            end
            [summaryLocal, metricsByStrategy] = run_prediction_sequence(subset, waferNames, signalType, cfg, resultDir);
            if ~isempty(summaryLocal)
                summaryLocal.recipe = repmat({recipeName}, height(summaryLocal), 1);
                summaryLocal.lot = repmat({lotName}, height(summaryLocal), 1);
                summaryLocal.parity = repmat({parity}, height(summaryLocal), 1);
                summaryLocal.signal = repmat({signalType}, height(summaryLocal), 1);
                summaryTbl = [summaryTbl; summaryLocal]; %#ok<AGROW>
            end
            save(fullfile(resultDir, 'metrics_by_strategy.mat'), 'metricsByStrategy');
        end
    end
end

if ~isempty(summaryTbl)
    writetable(summaryTbl, fullfile(outputDir, sprintf('summary_%s.csv', recipeName)));
end
end

function name = get_last_path_name(p)
[~, name] = fileparts(p);
end
