function filesSorted = collect_wafer_files(lotPath)
f1 = dir(fullfile(lotPath, '*.wmaftcor'));
f2 = dir(fullfile(lotPath, '*.log'));
allf = [f1; f2];
if isempty(allf)
    filesSorted = {};
    return;
end
[~, ord] = sort([allf.datenum], 'ascend');
allf = allf(ord);
filesSorted = cellfun(@(n) fullfile(lotPath, n), {allf.name}, 'UniformOutput', false);
end
