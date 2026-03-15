% ================================================================
% wafer_lot_consistency.m
% 硅片面型 lot 内一致性分析 —— 基于逐片预测误差
%
% 分析逻辑（模拟实际工况）:
%   按文件修改时间排序后分奇偶工位，每个工位独立分析。
%   对工位内第 X 片 (X>=2):
%     1) 前 X-1 片的平均面型 → 扣除其一阶拟合 → 残差均值 Y
%     2) 第 X 片在 3 个固定采样点 (r=75mm, 0/120/240 deg) 取 z 值
%        → 拟合一阶平面
%     3) 预测面型 = 该一阶平面 + Y
%     4) 误差面型 = 实际第X片 - 预测面型
%   收集所有误差面型，做一致性统计（Std/PV/径向等）。
%
% 批量处理：总文件夹 -> recipe子文件夹 -> lot子文件夹
% MATLAB R2020b 兼容
% ================================================================
clear; clc; close all;

%% ===== 0. 全局参数 =====
radius       = 150;      % 晶圆半径 (mm)
grid_res     = 1;        % 网格分辨率 (mm)
coef_xy      = 1e3;      % m -> mm
coef_z       = 1e9;      % m -> nm
sigma_thresh = 3;        % 跳点剔除阈值 (sigma)

% 3 个固定采样点: r=75mm, theta = 0, 120, 240 deg
sample_r     = 75;       % mm
sample_ang   = [0, 120, 240] * pi / 180;   % rad
sample_pts   = [sample_r * cos(sample_ang)', ...
                sample_r * sin(sample_ang)'];   % [3 x 2]

% 构建插值网格
[X_grid, Y_grid] = meshgrid(-radius:grid_res:radius, -radius:grid_res:radius);
mask = (X_grid.^2 + Y_grid.^2) <= radius^2;
n_rows = size(X_grid, 1);
n_cols = size(X_grid, 2);

% 网格的 x/y 向量（用于 interp2）
xvec = -radius:grid_res:radius;
yvec = -radius:grid_res:radius;

%% ===== 1. 选择总文件夹并遍历 =====
root_dir = uigetdir('', '选择总文件夹（包含recipe子文件夹）');
if isequal(root_dir, 0); error('未选择文件夹，程序终止。'); end

results_root = fullfile(root_dir, 'Consistency_Results');
if ~exist(results_root, 'dir'); mkdir(results_root); end

% 遍历 recipe 子文件夹
recipe_items = dir(root_dir);
recipe_items = recipe_items([recipe_items.isdir]);
recipe_items = recipe_items(~ismember({recipe_items.name}, {'.', '..', 'Consistency_Results'}));

total_lots = 0;

for ri = 1:length(recipe_items)
    recipe_name = recipe_items(ri).name;
    recipe_dir  = fullfile(root_dir, recipe_name);

    % 遍历 lot 子文件夹
    lot_items = dir(recipe_dir);
    lot_items = lot_items([lot_items.isdir]);
    lot_items = lot_items(~ismember({lot_items.name}, {'.', '..'}));

    for li = 1:length(lot_items)
        lot_name = lot_items(li).name;
        lot_dir  = fullfile(recipe_dir, lot_name);

        % 查找数据文件（.wmaftcor 和 .log）
        files_wma = dir(fullfile(lot_dir, '*.wmaftcor'));
        files_log = dir(fullfile(lot_dir, '*.log'));
        files_all = [files_wma; files_log];

        if isempty(files_all)
            fprintf('[跳过] 空文件夹: %s/%s\n', recipe_name, lot_name);
            continue;
        end

        % 按文件修改时间排序 → 奇偶工位依据
        [~, sort_idx] = sort([files_all.datenum]);
        files_all = files_all(sort_idx);
        num_files = length(files_all);

        fprintf('\n============================================================\n');
        fprintf('  处理: Recipe = %s,  Lot = %s  (%d 个文件)\n', ...
            recipe_name, lot_name, num_files);
        fprintf('============================================================\n');

        % 创建保存目录
        save_dir = fullfile(results_root, recipe_name, lot_name);
        if ~exist(save_dir, 'dir'); mkdir(save_dir); end

        % -----------------------------------------------------------
        %  读取与预处理所有晶圆到公共网格
        % -----------------------------------------------------------
        Z_all       = NaN(n_rows, n_cols, num_files);
        wafer_names = cell(num_files, 1);
        valid_flags = true(num_files, 1);

        for wi = 1:num_files
            filepath = fullfile(lot_dir, files_all(wi).name);
            wafer_names{wi} = files_all(wi).name;

            try
                raw = fun_read_wafer_map(filepath);
            catch ME
                fprintf('  [警告] 读取失败: %s  (%s)\n', files_all(wi).name, ME.message);
                valid_flags(wi) = false;
                continue;
            end

            x = raw(:, 6) * coef_xy;   % m -> mm
            y = raw(:, 7) * coef_xy;
            z = raw(:, 8) * coef_z;    % m -> nm

            % 3sigma 跳点剔除（二次曲面预拟合）
            A_pre = [x.^2, y.^2, x, y, ones(length(x), 1)];
            c_pre = A_pre \ z;
            res_pre = z - A_pre * c_pre;
            keep = abs(res_pre) <= sigma_thresh * std(res_pre);
            x = x(keep);  y = y(keep);  z = z(keep);

            if length(x) < 50
                fprintf('  [警告] 有效点不足: %s (%d 点)\n', files_all(wi).name, length(x));
                valid_flags(wi) = false;
                continue;
            end

            % 插值到公共网格
            F_interp = scatteredInterpolant(x, y, z, 'natural', 'none');
            Z_grid = F_interp(X_grid, Y_grid);
            Z_grid(~mask) = NaN;
            Z_all(:, :, wi) = Z_grid;

            fprintf('  [%2d/%2d] %s  -> OK\n', wi, num_files, files_all(wi).name);
        end

        % 移除无效晶圆
        invalid = ~valid_flags;
        if any(invalid)
            Z_all(:, :, invalid)  = [];
            wafer_names(invalid)  = [];
        end
        num_wafers = size(Z_all, 3);
        fprintf('  有效晶圆: %d / %d\n', num_wafers, num_files);

        % 奇偶分组（按排序后序号）
        idx_odd  = 1:2:num_wafers;
        idx_even = 2:2:num_wafers;

        % ===========================================================
        %  奇数工位分析
        % ===========================================================
        fprintf('\n  --- 奇数工位分析 (n=%d) ---\n', length(idx_odd));
        [err_odd, info_odd] = compute_prediction_errors( ...
            Z_all(:,:,idx_odd), X_grid, Y_grid, mask, ...
            xvec, yvec, sample_pts);

        if ~isempty(err_odd)
            stats_odd = run_error_consistency(err_odd, X_grid, Y_grid, mask, ...
                radius, '奇数工位', save_dir, 'F1_Odd', 'F2_Odd');
        else
            stats_odd = [];
            fprintf('  [跳过] 奇数工位有效误差面型不足。\n');
        end

        % ===========================================================
        %  偶数工位分析
        % ===========================================================
        fprintf('\n  --- 偶数工位分析 (n=%d) ---\n', length(idx_even));
        [err_even, info_even] = compute_prediction_errors( ...
            Z_all(:,:,idx_even), X_grid, Y_grid, mask, ...
            xvec, yvec, sample_pts);

        if ~isempty(err_even)
            stats_even = run_error_consistency(err_even, X_grid, Y_grid, mask, ...
                radius, '偶数工位', save_dir, 'F1_Even', 'F2_Even');
        else
            stats_even = [];
            fprintf('  [跳过] 偶数工位有效误差面型不足。\n');
        end

        % ===========================================================
        %  奇偶对比图
        % ===========================================================
        if ~isempty(stats_odd) && ~isempty(stats_even)
            plot_oddeven_compare(stats_odd, stats_even, ...
                X_grid, Y_grid, mask, save_dir);
        end

        % ===========================================================
        %  生成文字报告
        % ===========================================================
        write_report(save_dir, recipe_name, lot_name, ...
            num_wafers, num_files, wafer_names, ...
            idx_odd, idx_even, ...
            info_odd, info_even, stats_odd, stats_even);

        total_lots = total_lots + 1;
        fprintf('  >>> 完成: %s / %s\n', recipe_name, lot_name);
    end
end

fprintf('\n============================================================\n');
fprintf('  全部处理完毕。共处理 %d 个 lot。\n', total_lots);
fprintf('  结果保存路径: %s\n', results_root);
fprintf('============================================================\n');


%% =====================================================================
%%  局部函数
%% =====================================================================

% ----------------------------------------------------------------------
%  核心：计算工位内逐片预测误差面型
%  输入: Z_grp(:,:,N) — 该工位 N 片的网格面型
%  输出: Err_all(:,:,M) — 第2片起的误差面型 (M <= N-1)
%        info — 逐片详细记录 struct array
% ----------------------------------------------------------------------
function [Err_all, info] = compute_prediction_errors( ...
        Z_grp, X_grid, Y_grid, mask, xvec, yvec, sample_pts)

    N = size(Z_grp, 3);

    % 预分配 info 结构体数组
    empty_info = struct('idx', [], 'c1_sample', [], ...
                        'err_rms', [], 'err_pv', [], 'valid', []);
    info = repmat(empty_info, 0, 1);

    if N < 2
        Err_all = [];
        return;
    end

    Err_all = NaN(size(X_grid, 1), size(X_grid, 2), N - 1);
    valid_count = 0;

    for k = 2:N
        % 1) 前 k-1 片的平均面型
        if k == 2
            MeanSurf = Z_grp(:, :, 1);
        else
            MeanSurf = squeeze(mean(Z_grp(:, :, 1:k-1), 3, 'omitnan'));
        end

        % 2) 对平均面型做一阶拟合并扣除 → 残差均值 Y
        valid_mean = mask & ~isnan(MeanSurf);
        xm = X_grid(valid_mean);
        ym = Y_grid(valid_mean);
        zm = MeanSurf(valid_mean);

        if length(xm) < 10
            new_info.idx = k;
            new_info.valid = false;
            new_info.c1_sample = [NaN NaN NaN];
            new_info.err_rms = NaN;
            new_info.err_pv  = NaN;
            info = [info; new_info];  %#ok<AGROW>
            fprintf('    工位内第 %d 片: 平均面型有效点不足，跳过。\n', k);
            continue;
        end

        A_mean = [xm, ym, ones(length(xm), 1)];
        c_mean = A_mean \ zm;
        Fit1_mean = c_mean(1) .* X_grid + c_mean(2) .* Y_grid + c_mean(3);
        Fit1_mean(~mask) = NaN;
        Y_residual = MeanSurf - Fit1_mean;   % 残差均值

        % 3) 第 k 片：在 3 个采样点取 z 值
        Z_k = Z_grp(:, :, k);
        z_samples = NaN(3, 1);
        for sp = 1:3
            px = sample_pts(sp, 1);
            py = sample_pts(sp, 2);
            % interp2(X_vec, Y_vec, V, Xq, Yq)
            val = interp2(xvec, yvec, Z_k, px, py, 'linear');
            if isnan(val)
                val = interp2(xvec, yvec, Z_k, px, py, 'nearest');
            end
            z_samples(sp) = val;
        end

        if any(isnan(z_samples))
            new_info.idx = k;
            new_info.valid = false;
            new_info.c1_sample = [NaN NaN NaN];
            new_info.err_rms = NaN;
            new_info.err_pv  = NaN;
            info = [info; new_info];  %#ok<AGROW>
            fprintf('    工位内第 %d 片: 采样点存在 NaN，跳过。\n', k);
            continue;
        end

        % 4) 3 点拟合一阶平面: z = a*x + b*y + c  (精确解，3点3未知数)
        A_sample = [sample_pts, ones(3, 1)];
        c_sample = A_sample \ z_samples;

        Fit1_sample = c_sample(1) .* X_grid + c_sample(2) .* Y_grid + c_sample(3);
        Fit1_sample(~mask) = NaN;

        % 5) 预测面型 = 一阶平面 + Y_residual
        Predicted = Fit1_sample + Y_residual;

        % 6) 误差面型
        Error_k = Z_k - Predicted;

        valid_count = valid_count + 1;
        Err_all(:, :, valid_count) = Error_k;

        % 记录信息
        err_vals = Error_k(mask & ~isnan(Error_k));
        new_info.idx        = k;
        new_info.valid      = true;
        new_info.c1_sample  = c_sample';
        new_info.err_rms    = local_rms(err_vals);
        new_info.err_pv     = max(err_vals) - min(err_vals);
        info = [info; new_info];  %#ok<AGROW>

        fprintf('    工位内第 %d 片: 误差 RMS = %.2f nm,  PV = %.2f nm\n', ...
            k, new_info.err_rms, new_info.err_pv);
    end

    % 截断未使用的切片
    if valid_count > 0
        Err_all = Err_all(:, :, 1:valid_count);
    else
        Err_all = [];
    end
end


% ----------------------------------------------------------------------
%  误差面型一致性分析 + 可视化
% ----------------------------------------------------------------------
function stats = run_error_consistency(Err_all, X_grid, Y_grid, mask, ...
        radius, title_prefix, save_dir, tag_fig1, tag_fig2)

    num_e = size(Err_all, 3);

    % ---- 像素级批次统计（沿第3维），squeeze 兼容 R2020b ----
    Mean_map = squeeze(mean(Err_all, 3, 'omitnan'));
    Std_map  = squeeze(std(Err_all, 0, 3, 'omitnan'));
    Max_map  = squeeze(max(Err_all, [], 3, 'omitnan'));
    Min_map  = squeeze(min(Err_all, [], 3, 'omitnan'));
    PV_map   = Max_map - Min_map;

    % ---- 逐片指标 ----
    wafer_RMS  = zeros(num_e, 1);
    wafer_PV   = zeros(num_e, 1);
    wafer_dRMS = zeros(num_e, 1);     % 与误差均值的偏差 RMS

    for ei = 1:num_e
        surf_i = Err_all(:, :, ei);
        valid  = mask & ~isnan(surf_i);
        vals_i = surf_i(valid);
        wafer_RMS(ei) = local_rms(vals_i);
        wafer_PV(ei)  = max(vals_i) - min(vals_i);

        valid2 = valid & ~isnan(Mean_map);
        if sum(valid2(:)) > 10
            wafer_dRMS(ei) = local_rms(surf_i(valid2) - Mean_map(valid2));
        end
    end

    % ---- 汇总统计 ----
    std_vals = Std_map(mask & ~isnan(Std_map));
    pv_vals  = PV_map(mask & ~isnan(PV_map));

    stats.Mean_map   = Mean_map;
    stats.Std_map    = Std_map;
    stats.PV_map     = PV_map;
    stats.wafer_RMS  = wafer_RMS;
    stats.wafer_PV   = wafer_PV;
    stats.wafer_dRMS = wafer_dRMS;
    stats.num_errors = num_e;
    stats.Std_mean   = mean(std_vals);
    stats.Std_max    = max(std_vals);
    stats.Std_median = median(std_vals);
    stats.PV_mean    = mean(pv_vals);
    stats.PV_max     = max(pv_vals);

    % ---- 径向一致性剖面 ----
    R_map = sqrt(X_grid.^2 + Y_grid.^2);
    r_edges = 0:5:radius;
    r_centers = r_edges(1:end-1) + 2.5;
    n_bins = length(r_centers);
    radial_std_mean = zeros(n_bins, 1);
    radial_std_std  = zeros(n_bins, 1);

    for bi = 1:n_bins
        ring = mask & (R_map >= r_edges(bi)) & (R_map < r_edges(bi+1));
        ring_std = Std_map(ring & ~isnan(Std_map));
        if ~isempty(ring_std)
            radial_std_mean(bi) = mean(ring_std);
            radial_std_std(bi)  = std(ring_std);
        end
    end
    stats.r_centers       = r_centers;
    stats.radial_std_mean = radial_std_mean;
    stats.radial_std_std  = radial_std_std;

    % ============== 绘图 ==============

    % ---- Figure 1: 空间统计图 (2x2) ----
    fig1 = figure('Position', [50 50 1300 1000], 'Visible', 'off');

    subplot(2, 2, 1);
    pcolor(X_grid, Y_grid, Mean_map); shading interp; colorbar;
    axis equal tight; set(gca, 'FontSize', 10);
    title(sprintf('%s - 平均误差面型 (nm)', title_prefix), 'FontSize', 12);
    xlabel('X (mm)'); ylabel('Y (mm)');

    subplot(2, 2, 2);
    pcolor(X_grid, Y_grid, Std_map); shading interp; colorbar;
    axis equal tight; set(gca, 'FontSize', 10);
    title(sprintf('%s - 误差Std分布 (nm)', title_prefix), 'FontSize', 12);
    xlabel('X (mm)'); ylabel('Y (mm)');

    subplot(2, 2, 3);
    pcolor(X_grid, Y_grid, PV_map); shading interp; colorbar;
    axis equal tight; set(gca, 'FontSize', 10);
    title(sprintf('%s - 误差PV分布 (nm)', title_prefix), 'FontSize', 12);
    xlabel('X (mm)'); ylabel('Y (mm)');

    subplot(2, 2, 4);
    histogram(std_vals, 50, 'FaceColor', [0.3 0.6 0.9], 'EdgeColor', 'w');
    set(gca, 'FontSize', 10);
    xlabel('Std (nm)'); ylabel('像素数');
    title(sprintf('%s - Std柱状图', title_prefix), 'FontSize', 12);
    xl = xlim;
    text(xl(2)*0.95, max(ylim)*0.90, ...
        sprintf('Mean = %.2f nm\nMax = %.2f nm\nMedian = %.2f nm', ...
        stats.Std_mean, stats.Std_max, stats.Std_median), ...
        'HorizontalAlignment', 'right', 'FontSize', 9, ...
        'BackgroundColor', 'w', 'EdgeColor', [0.5 0.5 0.5]);

    sgtitle(sprintf('%s — 预测误差空间一致性', title_prefix), ...
        'FontSize', 14, 'FontWeight', 'bold');
    saveas(fig1, fullfile(save_dir, sprintf('%s_SpatialStats.png', tag_fig1)));
    close(fig1);

    % ---- Figure 2: 逐片指标 + 径向剖面 (2x2) ----
    fig2 = figure('Position', [50 50 1300 1000], 'Visible', 'off');
    seq = 1:num_e;

    subplot(2, 2, 1);
    bar(seq, wafer_RMS, 'FaceColor', [0.2 0.6 0.8]);
    set(gca, 'FontSize', 10);
    xlabel('序号（工位内，从第2片起）'); ylabel('误差 RMS (nm)');
    title(sprintf('%s - 逐片误差 RMS', title_prefix), 'FontSize', 12);
    hold on;
    yline(mean(wafer_RMS), 'r--', sprintf('Mean=%.1f', mean(wafer_RMS)), ...
        'LineWidth', 1.2, 'FontSize', 9, 'LabelHorizontalAlignment', 'left');
    hold off;

    subplot(2, 2, 2);
    bar(seq, wafer_PV, 'FaceColor', [0.9 0.5 0.2]);
    set(gca, 'FontSize', 10);
    xlabel('序号（工位内，从第2片起）'); ylabel('误差 PV (nm)');
    title(sprintf('%s - 逐片误差 PV', title_prefix), 'FontSize', 12);
    hold on;
    yline(mean(wafer_PV), 'r--', sprintf('Mean=%.1f', mean(wafer_PV)), ...
        'LineWidth', 1.2, 'FontSize', 9, 'LabelHorizontalAlignment', 'left');
    hold off;

    subplot(2, 2, 3);
    plot(seq, wafer_RMS, 'b-o', 'LineWidth', 1.2, 'MarkerFaceColor', [0.2 0.6 0.8]);
    set(gca, 'FontSize', 10);
    xlabel('序号（工位内，从第2片起）'); ylabel('误差 RMS (nm)');
    title(sprintf('%s - 误差RMS趋势', title_prefix), 'FontSize', 12);
    grid on;

    subplot(2, 2, 4);
    fill([r_centers, fliplr(r_centers)], ...
         [radial_std_mean' + radial_std_std', ...
          fliplr(radial_std_mean' - radial_std_std')], ...
         [0.8 0.9 1.0], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    hold on;
    plot(r_centers, radial_std_mean, 'b-', 'LineWidth', 1.5);
    hold off;
    set(gca, 'FontSize', 10);
    xlabel('径向距离 (mm)'); ylabel('误差 Std (nm)');
    title(sprintf('%s - 径向一致性剖面 (Mean +/- Std)', title_prefix), 'FontSize', 12);
    xlim([0 radius]);

    sgtitle(sprintf('%s — 逐片误差指标与径向剖面', title_prefix), ...
        'FontSize', 14, 'FontWeight', 'bold');
    saveas(fig2, fullfile(save_dir, sprintf('%s_PerWafer.png', tag_fig2)));
    close(fig2);
end


% ----------------------------------------------------------------------
%  奇偶对比图
% ----------------------------------------------------------------------
function plot_oddeven_compare(stats_odd, stats_even, ...
        X_grid, Y_grid, mask, save_dir)

    fig = figure('Position', [50 50 1500 900], 'Visible', 'off');

    % Row 1: Std 空间分布对比
    subplot(2, 3, 1);
    pcolor(X_grid, Y_grid, stats_odd.Std_map); shading interp; colorbar;
    axis equal tight; title(sprintf('奇数-误差Std (Mean=%.2f nm)', stats_odd.Std_mean), 'FontSize', 11);

    subplot(2, 3, 2);
    pcolor(X_grid, Y_grid, stats_even.Std_map); shading interp; colorbar;
    axis equal tight; title(sprintf('偶数-误差Std (Mean=%.2f nm)', stats_even.Std_mean), 'FontSize', 11);

    subplot(2, 3, 3);
    Diff_std = stats_odd.Std_map - stats_even.Std_map;
    pcolor(X_grid, Y_grid, Diff_std); shading interp; colorbar;
    axis equal tight;
    diff_vals = Diff_std(mask & ~isnan(Diff_std));
    title(sprintf('Std差异(奇-偶) RMS=%.2f nm', local_rms(diff_vals)), 'FontSize', 11);

    % Row 2: 逐片 RMS 对比 + 径向剖面对比 + 汇总柱状图
    subplot(2, 3, 4);
    n_o = stats_odd.num_errors;
    n_e = stats_even.num_errors;
    hold on;
    plot(1:n_o, stats_odd.wafer_RMS,  'b-o', 'LineWidth', 1.2, ...
        'MarkerFaceColor', [0.3 0.6 0.9]);
    plot(1:n_e, stats_even.wafer_RMS, 'r-s', 'LineWidth', 1.2, ...
        'MarkerFaceColor', [0.9 0.5 0.2]);
    hold off;
    xlabel('序号（工位内，从第2片起）'); ylabel('误差 RMS (nm)');
    title('逐片误差RMS趋势', 'FontSize', 11);
    legend(sprintf('奇数(n=%d)', n_o), sprintf('偶数(n=%d)', n_e), 'Location', 'best');
    set(gca, 'FontSize', 10); grid on;

    subplot(2, 3, 5);
    hold on;
    plot(stats_odd.r_centers,  stats_odd.radial_std_mean,  'b-', 'LineWidth', 1.5);
    plot(stats_even.r_centers, stats_even.radial_std_mean, 'r-', 'LineWidth', 1.5);
    hold off;
    xlabel('径向距离 (mm)'); ylabel('误差 Std (nm)');
    title('径向一致性剖面对比', 'FontSize', 11);
    legend('奇数', '偶数', 'Location', 'best');
    set(gca, 'FontSize', 10); xlim([0 max(stats_odd.r_centers)]);

    subplot(2, 3, 6);
    bar_data = [stats_odd.Std_mean,  stats_even.Std_mean;
                stats_odd.Std_max,   stats_even.Std_max;
                mean(stats_odd.wafer_RMS), mean(stats_even.wafer_RMS);
                mean(stats_odd.wafer_PV),  mean(stats_even.wafer_PV)];
    b = bar(bar_data);
    b(1).FaceColor = [0.3 0.6 0.9]; b(2).FaceColor = [0.9 0.5 0.2];
    set(gca, 'XTickLabel', {'Mean(Std)', 'Max(Std)', 'Mean(RMS)', 'Mean(PV)'}, ...
        'FontSize', 9);
    ylabel('nm'); title('汇总指标对比', 'FontSize', 11);
    legend('奇数', '偶数', 'Location', 'northwest');

    sgtitle('奇偶工位预测误差对比', 'FontSize', 14, 'FontWeight', 'bold');
    saveas(fig, fullfile(save_dir, 'F3_OddEven_Compare.png'));
    close(fig);
end


% ----------------------------------------------------------------------
%  生成纯文本报告
% ----------------------------------------------------------------------
function write_report(save_dir, recipe_name, lot_name, ...
        num_wafers, num_files, wafer_names, ...
        idx_odd, idx_even, ...
        info_odd, info_even, stats_odd, stats_even)

    fpath = fullfile(save_dir, 'Consistency_Report.txt');
    fid = fopen(fpath, 'w', 'n', 'UTF-8');
    pw = @(fmt, varargin) fprintf(fid, [fmt '\n'], varargin{:});

    pw('================================================================');
    pw('  硅片面型 LOT 内一致性分析报告（逐片预测误差）');
    pw('================================================================');
    pw('');
    pw('生成时间:  %s', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    pw('Recipe:    %s', recipe_name);
    pw('Lot:       %s', lot_name);
    pw('文件总数:  %d', num_files);
    pw('有效晶圆:  %d', num_wafers);
    pw('奇数工位:  %d 片 (全局序号 %s)', length(idx_odd),  mat2str(idx_odd));
    pw('偶数工位:  %d 片 (全局序号 %s)', length(idx_even), mat2str(idx_even));
    pw('');
    pw('采样点:  r = 75 mm,  角度 = 0, 120, 240 deg');
    pw('');
    pw('分析方法:');
    pw('  对工位内第 X 片 (X>=2):');
    pw('    1) 前 X-1 片平均面型 -> 扣除一阶 -> 残差均值 Y');
    pw('    2) 第 X 片在 3 采样点取值 -> 拟合一阶平面');
    pw('    3) 预测 = 一阶平面 + Y;  误差 = 实际 - 预测');
    pw('');

    % ---- 晶圆文件清单 ----
    pw('================================================================');
    pw('  晶圆文件清单');
    pw('================================================================');
    pw('  全局序号  工位  工位内序号  文件名');
    pw('  --------  ----  ----------  ------------------------------------');
    for wi = 1:num_wafers
        if mod(wi, 2) == 1
            station = '奇';
            grp_idx = ceil(wi / 2);
        else
            station = '偶';
            grp_idx = wi / 2;
        end
        name_str = wafer_names{wi};
        if length(name_str) > 36; name_str = [name_str(1:33) '...']; end
        pw('  %8d  %s    %10d  %s', wi, station, grp_idx, name_str);
    end
    pw('');

    % ---- 奇数工位 ----
    pw('================================================================');
    pw('  Part 1: 奇数工位 — 逐片预测误差');
    pw('================================================================');
    write_station_section(fid, pw, info_odd, stats_odd, '奇数', idx_odd, wafer_names);

    % ---- 偶数工位 ----
    pw('================================================================');
    pw('  Part 2: 偶数工位 — 逐片预测误差');
    pw('================================================================');
    write_station_section(fid, pw, info_even, stats_even, '偶数', idx_even, wafer_names);

    % ---- 奇偶对比 ----
    pw('================================================================');
    pw('  Part 3: 奇偶工位对比');
    pw('================================================================');
    if ~isempty(stats_odd) && ~isempty(stats_even)
        pw('');
        pw('  指标                  奇数工位       偶数工位');
        pw('  --------------------  -------------  -------------');
        pw('  误差Std均值 (nm)      %13.2f  %13.2f', stats_odd.Std_mean,  stats_even.Std_mean);
        pw('  误差Std最大 (nm)      %13.2f  %13.2f', stats_odd.Std_max,   stats_even.Std_max);
        pw('  误差PV均值  (nm)      %13.2f  %13.2f', stats_odd.PV_mean,   stats_even.PV_mean);
        pw('  误差PV最大  (nm)      %13.2f  %13.2f', stats_odd.PV_max,    stats_even.PV_max);
        pw('  逐片RMS均值 (nm)      %13.2f  %13.2f', mean(stats_odd.wafer_RMS), mean(stats_even.wafer_RMS));
        pw('  逐片PV均值  (nm)      %13.2f  %13.2f', mean(stats_odd.wafer_PV),  mean(stats_even.wafer_PV));
        pw('  误差片数               %13d  %13d', stats_odd.num_errors, stats_even.num_errors);
    else
        pw('  数据不足，无法对比。');
    end
    pw('');

    % ---- 输出文件列表 ----
    pw('================================================================');
    pw('  输出文件列表');
    pw('================================================================');
    pw('  F1_Odd_SpatialStats.png    奇数工位 - 误差空间统计');
    pw('  F2_Odd_PerWafer.png        奇数工位 - 逐片误差指标');
    pw('  F1_Even_SpatialStats.png   偶数工位 - 误差空间统计');
    pw('  F2_Even_PerWafer.png       偶数工位 - 逐片误差指标');
    pw('  F3_OddEven_Compare.png     奇偶对比');
    pw('  Consistency_Report.txt     本报告');
    pw('');
    pw('========================== 报告结束 ==========================');

    fclose(fid);
    fprintf('  报告已保存: %s\n', fpath);
end


% 写入单工位的详细报告段
function write_station_section(fid, pw, info, stats, station_name, ...
        global_idx, wafer_names)
    pw('');
    if isempty(info)
        pw('  该工位晶圆不足 2 片，无误差数据。');
        pw('');
        return;
    end

    pw('  (工位内第 1 片为基准，不产生误差)');
    pw('');
    pw('  工位内序号  全局序号  文件名                                误差RMS(nm)  误差PV(nm)  状态');
    pw('  ----------  --------  ------------------------------------  -----------  ----------  ------');
    for ii = 1:length(info)
        grp_seq  = info(ii).idx;              % 工位内序号
        glob_seq = global_idx(grp_seq);       % 全局序号
        name_str = wafer_names{glob_seq};
        if length(name_str) > 36; name_str = [name_str(1:33) '...']; end
        if info(ii).valid
            pw('  %10d  %8d  %-36s  %11.2f  %10.2f  有效', ...
                grp_seq, glob_seq, name_str, info(ii).err_rms, info(ii).err_pv);
        else
            pw('  %10d  %8d  %-36s  %11s  %10s  跳过', ...
                grp_seq, glob_seq, name_str, '-', '-');
        end
    end
    pw('');

    if ~isempty(stats)
        pw('  汇总统计 (基于 %d 个有效误差面型):', stats.num_errors);
        pw('    像素级Std  -- 均值: %.2f nm,  最大: %.2f nm,  中位数: %.2f nm', ...
            stats.Std_mean, stats.Std_max, stats.Std_median);
        pw('    像素级PV   -- 均值: %.2f nm,  最大: %.2f nm', stats.PV_mean, stats.PV_max);
        pw('    逐片RMS    (Mean +/- Std): %.2f +/- %.2f nm,  范围: [%.2f, %.2f]', ...
            mean(stats.wafer_RMS), std(stats.wafer_RMS), min(stats.wafer_RMS), max(stats.wafer_RMS));
        pw('    逐片PV     (Mean +/- Std): %.2f +/- %.2f nm,  范围: [%.2f, %.2f]', ...
            mean(stats.wafer_PV), std(stats.wafer_PV), min(stats.wafer_PV), max(stats.wafer_PV));
    else
        pw('  有效误差面型不足，未生成统计。');
    end
    pw('');
end


% ----------------------------------------------------------------------
%  本地 RMS 函数（避免依赖 Signal Processing Toolbox）
% ----------------------------------------------------------------------
function val = local_rms(x)
    val = sqrt(mean(x(:).^2));
end
