Wafer edge-field prediction benchmark (MATLAB R2020b)

Entry point:
  run_edge_field_pipeline(rootDir, outputDir)

Data assumption:
  rootDir/recipe_folder/lot_folder/*.wmaftcor or *.log
  Each file is read by: data = fun_read_wafer_map(filepath)
  data(:,6)=x, data(:,7)=y, data(:,8)=z (scatter points)

Main behaviors:
  - Sort wafers by file modified time.
  - Odd/even wafer separation by sorted index.
  - Edge region: 100 < r <= 150 mm.
  - Notch exclusion: default center=270 deg, half-width=10 deg.
  - Outlier handling: >3 sigma replaced by scattered interpolation.
  - First-order removal Y = X - Zplane using 3-point leveling:
      angles [0,120,240] deg at r=75 mm.

Strategies:
  1) mean_previous
  2) rbf_gaussian
  3) bspline_tensor
  4) zernike
  5) fourier_radial (additional candidate)

Outputs:
  - Per recipe/lot/parity/signal summary.csv
  - Per strategy: metrics CSV + PNG figures
  - Aggregate summary_all_recipes.csv
  - workspace_all.mat
