Wafer edge-field prediction benchmark (MATLAB R2020b)

Quick start (recommended)
-------------------------
1) Ensure `fun_read_wafer_map.m` is on MATLAB path.
2) Open `run_edge_field_pipeline_main.m`.
3) Set:
   ROOT_DIR   = your wafer root folder
   OUTPUT_DIR = folder for results
4) Run `run_edge_field_pipeline_main.m`.

Alternative entry
-----------------
Call directly in MATLAB command window:
  run_edge_field_pipeline(rootDir, outputDir)

Data assumption
---------------
rootDir/recipe_folder/lot_folder/*.wmaftcor or *.log
Each file is read by: data = fun_read_wafer_map(filepath)
data(:,6)=x, data(:,7)=y, data(:,8)=z (scatter points)

Main behaviors
--------------
- Sort wafers by file modified time.
- Odd/even wafer separation by sorted index.
- Edge region: 100 < r <= 150 mm.
- Notch exclusion: default center=270 deg, half-width=10 deg.
- Outlier handling: >3 sigma replaced by scattered interpolation.
- First-order removal Y = X - Zplane using 3-point leveling:
    angles [0,120,240] deg at r=75 mm.
- Rolling prediction rule: for wafer X (>1), train with wafers 1..X-1 then predict X.

Strategies
----------
1) mean_previous
2) rbf_gaussian
3) bspline_tensor
4) zernike
5) fourier_radial (additional candidate)

How to know which strategy is best
----------------------------------
After run completes, check:
- `best_strategy_overall.csv`: best strategy in each (recipe/lot/parity/signal) group.
- `strategy_rank_overall.csv`: global ranking by mean RMSE (lower is better).
- `strategy_report.txt`: text report including top strategy and ranking details.

Outputs
-------
- Per recipe/lot/parity/signal summary.csv
- Per strategy: metrics CSV + PNG figures
- Average shape CSV/PNG
- Aggregate summary_all_recipes.csv
- best_strategy_overall.csv
- strategy_rank_overall.csv
- strategy_report.txt
- workspace_all.mat
