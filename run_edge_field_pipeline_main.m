% MATLAB R2020b one-click entry script.
% Usage:
%   1) Put this folder and fun_read_wafer_map.m in MATLAB path.
%   2) Edit ROOT_DIR/OUTPUT_DIR below.
%   3) Run this script.

ROOT_DIR = 'D:\\wafer_root';   % TODO: change to your root folder
OUTPUT_DIR = fullfile(ROOT_DIR, 'edge_field_results');

run_edge_field_pipeline(ROOT_DIR, OUTPUT_DIR);

% Print best strategies from generated summary:
summaryFile = fullfile(OUTPUT_DIR, 'best_strategy_overall.csv');
if exist(summaryFile, 'file')
    T = readtable(summaryFile);
    disp('==== Best strategy by group (recipe/lot/parity/signal) ====');
    disp(T(:, {'recipe','lot','parity','signal','best_strategy','best_rmse'}));
else
    warning('best_strategy_overall.csv not found.');
end
