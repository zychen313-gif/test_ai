function run_edge_field_pipeline(rootDir, outputDir)
% MATLAB R2020b
% Batch wafer edge-field prediction benchmark pipeline.

if nargin < 1 || isempty(rootDir)
    rootDir = pwd;
end
if nargin < 2 || isempty(outputDir)
    outputDir = fullfile(rootDir, 'edge_field_results');
end
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

cfg = default_config();
recipes = dir(rootDir);
recipes = recipes([recipes.isdir]);
recipes = recipes(~ismember({recipes.name}, {'.','..', 'edge_field_results'}));

allSummary = table();
for i = 1:numel(recipes)
    recipePath = fullfile(rootDir, recipes(i).name);
    fprintf('Processing recipe: %s\n', recipes(i).name);
    summaryTbl = process_recipe(recipePath, outputDir, cfg);
    if ~isempty(summaryTbl)
        allSummary = [allSummary; summaryTbl]; %#ok<AGROW>
    end
end

if ~isempty(allSummary)
    writetable(allSummary, fullfile(outputDir, 'summary_all_recipes.csv'));
end
save(fullfile(outputDir, 'workspace_all.mat'));
fprintf('Done. Results saved to: %s\n', outputDir);
end
