clc;
clear;
close all;

%% Runtime options
if_export = false;  % true: save PNG/EPS/PDF; false: display only

%% Geometrical-angle correction
% The raw angular calibration uses the discharge-off beam as its zero.
% Remove the wavelength-independent turning angle theta_0 = S/R so that
% the reported theta is measured relative to the corresponding channel axis.
S_mm = 2;
R20_mm = 20;
R40_mm = 40;
theta0_R20_deg = rad2deg(S_mm / R20_mm);
theta0_R40_deg = rad2deg(S_mm / R40_mm);

%% Paths
script_dir = fileparts(mfilename('fullpath'));
figs_root = fileparts(fileparts(script_dir));
data_dir = fullfile(figs_root, 'data', 'Fig2');
output_dir = fullfile(figs_root, 'Figs', 'Fig2');
if if_export && ~isfolder(output_dir)
    mkdir(output_dir);
end

%% Load data
stats = load(fullfile(data_dir, 'angular_dispersion_statistics.mat'));
stats_corrected = subtract_channel_turning_angle(stats, [1, 2], theta0_R20_deg);
stats_corrected = subtract_channel_turning_angle(stats_corrected, [3, 4], theta0_R40_deg);

%% Figure style
fs = 20;
fg = figure(2);
set(fg, 'Color', 'w', 'Units', 'pixels', ...
    'Position', [300, 120, 700, 560]);
ax_d = axes('Parent', fg);
hold(ax_d, 'on');

% Match the present Fig2_prx.m panel-(d) style.
if_errorbar = false;
if_shading = true;
if_connection = true;
shading_cfg.alpha = 0.2;
shading_cfg.expansion = 1.08;
shading_cfg.smoothing_iterations = 5;
psi10_size = 5;
psi20_size = 2;
reference_size = 2;
marker_cfg.psi20.marker = 'o';
marker_cfg.psi20.size = psi20_size;
marker_cfg.psi20.filled = true;
marker_cfg.psi10.marker = 'x';
marker_cfg.psi10.size = psi10_size;
marker_cfg.psi10.filled = false;
marker_cfg.reference.marker = 'o';
marker_cfg.reference.size = reference_size;
marker_cfg.reference.filled = true;

if if_shading
    plot_group_shading(ax_d, stats_corrected, [3, 4], [0.20, 0.35, 1.00], shading_cfg);
    plot_group_shading(ax_d, stats_corrected, [1, 2], [1.00, 0.20, 0.20], shading_cfg);
    plot_group_shading(ax_d, stats_corrected, 5, [0.25, 0.25, 0.25], shading_cfg);
end

if if_connection
    plot_centroid_connection(ax_d, stats_corrected, [3, 4], 'b');
    plot_centroid_connection(ax_d, stats_corrected, [1, 2], 'r');
    plot_centroid_connection(ax_d, stats_corrected, 5, 'k');
end

if if_errorbar
    plot_pooled_errorbars(ax_d, stats_corrected, [3, 4], 'b');
    plot_pooled_errorbars(ax_d, stats_corrected, [1, 2], 'r');
    plot_pooled_errorbars(ax_d, stats_corrected, 5, 'k');
end

plot_pair(ax_d, stats_corrected.lam1_5, stats_corrected.deg1_5, ...
    stats_corrected.lam2_5, stats_corrected.deg2_5, marker_cfg.reference, 'k');
plot_pair(ax_d, stats_corrected.lam1_3, stats_corrected.deg1_3, ...
    stats_corrected.lam2_3, stats_corrected.deg2_3, marker_cfg.psi20, 'b');
plot_pair(ax_d, stats_corrected.lam1_4, stats_corrected.deg1_4, ...
    stats_corrected.lam2_4, stats_corrected.deg2_4, marker_cfg.psi10, 'b');
plot_pair(ax_d, stats_corrected.lam1_1, stats_corrected.deg1_1, ...
    stats_corrected.lam2_1, stats_corrected.deg2_1, marker_cfg.psi20, 'r');
plot_pair(ax_d, stats_corrected.lam1_2, stats_corrected.deg1_2, ...
    stats_corrected.lam2_2, stats_corrected.deg2_2, marker_cfg.psi10, 'r');

xlim(ax_d, [0.5, 5]);
ylim(ax_d, [-6.0, 0.5]);
xlabel(ax_d, '$\lambda~(\mu\mathrm{m})$', 'Interpreter', 'latex');
ylabel(ax_d, '$\theta~(^{\circ})$', 'Interpreter', 'latex');
title(ax_d, 'Angular spectra');
set(ax_d, 'FontName', 'Times New Roman', 'FontSize', fs, ...
    'Box', 'on', 'TickDir', 'out');
add_grouped_legend(ax_d, fs, marker_cfg);

%% Statistics
group_table = build_group_statistics_table(stats_corrected);
pooled_table = build_pooled_statistics_table(stats_corrected);

fprintf('\nPer-group wavelength and angle statistics\n');
disp(group_table);

fprintf('\nPooled wavelength and angle statistics\n');
disp(pooled_table);

%% Export
if if_export
    output_png = fullfile(output_dir, 'Fig2_wavelength_energy_results_show.png');
    output_eps = fullfile(output_dir, 'Fig2_wavelength_energy_results_show.eps');
    output_pdf = fullfile(output_dir, 'Fig2_wavelength_energy_results_show.pdf');
    exportgraphics(fg, output_png, 'Resolution', 300);
    print(fg, output_eps, '-depsc', '-painters');
    exportgraphics(fg, output_pdf, 'ContentType', 'vector');
    fprintf('\nSaved figure to:\n%s\n%s\n%s\n', output_png, output_eps, output_pdf);
end

function stats = subtract_channel_turning_angle(stats, group_indices, theta0_deg)
    for group_idx = group_indices
        for spectral_group = 1:2
            field_name = sprintf('deg%d_%d', spectral_group, group_idx);
            stats.(field_name) = stats.(field_name) - theta0_deg;
        end
    end
end

function tbl = build_group_statistics_table(stats)
    group_defs = { ...
        1, 'R=20 mm', '20 psi'; ...
        2, 'R=20 mm', '10 psi'; ...
        3, 'R=40 mm', '20 psi'; ...
        4, 'R=40 mm', '10 psi'; ...
        5, 'No channel', 'No discharge'};

    n_groups = size(group_defs, 1);
    structure = strings(n_groups, 1);
    pressure = strings(n_groups, 1);
    nir_lambda_mean = nan(n_groups, 1);
    nir_lambda_std = nan(n_groups, 1);
    nir_lambda_n = nan(n_groups, 1);
    nir_theta_mean = nan(n_groups, 1);
    nir_theta_std = nan(n_groups, 1);
    nir_theta_n = nan(n_groups, 1);
    mir_lambda_mean = nan(n_groups, 1);
    mir_lambda_std = nan(n_groups, 1);
    mir_lambda_n = nan(n_groups, 1);
    mir_theta_mean = nan(n_groups, 1);
    mir_theta_std = nan(n_groups, 1);
    mir_theta_n = nan(n_groups, 1);

    for idx = 1:n_groups
        group_idx = group_defs{idx, 1};
        structure(idx) = string(group_defs{idx, 2});
        pressure(idx) = string(group_defs{idx, 3});

        nir_lambda = stats.(sprintf('lam1_%d', group_idx));
        nir_theta = stats.(sprintf('deg1_%d', group_idx));
        mir_lambda = stats.(sprintf('lam2_%d', group_idx));
        mir_theta = stats.(sprintf('deg2_%d', group_idx));

        [nir_lambda_mean(idx), nir_lambda_std(idx), nir_lambda_n(idx)] = summarize_vector(nir_lambda);
        [nir_theta_mean(idx), nir_theta_std(idx), nir_theta_n(idx)] = summarize_vector(nir_theta);
        [mir_lambda_mean(idx), mir_lambda_std(idx), mir_lambda_n(idx)] = summarize_vector(mir_lambda);
        [mir_theta_mean(idx), mir_theta_std(idx), mir_theta_n(idx)] = summarize_vector(mir_theta);
    end

    tbl = table(structure, pressure, ...
        nir_lambda_mean, nir_lambda_std, nir_lambda_n, ...
        nir_theta_mean, nir_theta_std, nir_theta_n, ...
        mir_lambda_mean, mir_lambda_std, mir_lambda_n, ...
        mir_theta_mean, mir_theta_std, mir_theta_n);
end

function tbl = build_pooled_statistics_table(stats)
    case_defs = { ...
        'R=20 mm', [1, 2]; ...
        'R=40 mm', [3, 4]; ...
        'No channel', 5};

    n_cases = size(case_defs, 1);
    structure = strings(n_cases, 1);
    nir_lambda_mean = nan(n_cases, 1);
    nir_lambda_std = nan(n_cases, 1);
    nir_lambda_n = nan(n_cases, 1);
    nir_theta_mean = nan(n_cases, 1);
    nir_theta_std = nan(n_cases, 1);
    nir_theta_n = nan(n_cases, 1);
    mir_lambda_mean = nan(n_cases, 1);
    mir_lambda_std = nan(n_cases, 1);
    mir_lambda_n = nan(n_cases, 1);
    mir_theta_mean = nan(n_cases, 1);
    mir_theta_std = nan(n_cases, 1);
    mir_theta_n = nan(n_cases, 1);
    delta_lambda = nan(n_cases, 1);
    delta_lambda_std = nan(n_cases, 1);
    delta_theta = nan(n_cases, 1);
    delta_theta_std = nan(n_cases, 1);
    angular_dispersion = nan(n_cases, 1);
    angular_dispersion_std = nan(n_cases, 1);

    for idx = 1:n_cases
        structure(idx) = string(case_defs{idx, 1});
        group_indices = case_defs{idx, 2};

        nir_lambda = collect_group_values(stats, group_indices, 'lam1');
        nir_theta = collect_group_values(stats, group_indices, 'deg1');
        mir_lambda = collect_group_values(stats, group_indices, 'lam2');
        mir_theta = collect_group_values(stats, group_indices, 'deg2');

        [nir_lambda_mean(idx), nir_lambda_std(idx), nir_lambda_n(idx)] = summarize_vector(nir_lambda);
        [nir_theta_mean(idx), nir_theta_std(idx), nir_theta_n(idx)] = summarize_vector(nir_theta);
        [mir_lambda_mean(idx), mir_lambda_std(idx), mir_lambda_n(idx)] = summarize_vector(mir_lambda);
        [mir_theta_mean(idx), mir_theta_std(idx), mir_theta_n(idx)] = summarize_vector(mir_theta);

        delta_lambda(idx) = mir_lambda_mean(idx) - nir_lambda_mean(idx);
        delta_lambda_std(idx) = hypot(mir_lambda_std(idx), nir_lambda_std(idx));
        delta_theta(idx) = mir_theta_mean(idx) - nir_theta_mean(idx);
        delta_theta_std(idx) = hypot(mir_theta_std(idx), nir_theta_std(idx));
        angular_dispersion(idx) = delta_theta(idx) / delta_lambda(idx);
        angular_dispersion_std(idx) = hypot( ...
            delta_theta_std(idx) / delta_lambda(idx), ...
            delta_theta(idx) * delta_lambda_std(idx) / delta_lambda(idx)^2);
    end

    tbl = table(structure, ...
        nir_lambda_mean, nir_lambda_std, nir_lambda_n, ...
        nir_theta_mean, nir_theta_std, nir_theta_n, ...
        mir_lambda_mean, mir_lambda_std, mir_lambda_n, ...
        mir_theta_mean, mir_theta_std, mir_theta_n, ...
        delta_lambda, delta_lambda_std, delta_theta, delta_theta_std, ...
        angular_dispersion, angular_dispersion_std);
end

function values = collect_group_values(stats, group_indices, field_prefix)
    values = [];
    for group_idx = group_indices
        values = [values, stats.(sprintf('%s_%d', field_prefix, group_idx))]; %#ok<AGROW>
    end
    values = values(isfinite(values));
end

function [value_mean, value_std, value_n] = summarize_vector(values)
    values = values(isfinite(values));
    value_n = numel(values);
    if value_n == 0
        value_mean = NaN;
        value_std = NaN;
    else
        value_mean = mean(values);
        value_std = std(values);
    end
end

function plot_pair(ax, lambda_1, theta_1, lambda_2, theta_2, style, color)
    marker_face = 'none';
    if style.filled
        marker_face = color;
    end
    plot(ax, lambda_1, theta_1, style.marker, 'Color', color, ...
        'MarkerFaceColor', marker_face, ...
        'MarkerSize', style.size, 'LineWidth', 1.3, ...
        'LineStyle', 'none');
    plot(ax, lambda_2, theta_2, style.marker, 'Color', color, ...
        'MarkerFaceColor', marker_face, ...
        'MarkerSize', style.size, 'LineWidth', 1.3, ...
        'LineStyle', 'none');
end

function add_grouped_legend(ax, fs, marker_cfg)
    x_move = 2.4;
    y_move = -3.3;
    rectangle(ax, 'Position', [0.58 + x_move, 0.76 + y_move, 1.85, 1.8], ...
        'FaceColor', 'w', 'EdgeColor', 'k', 'LineWidth', 0.8);
    x_marker_1 = 0.72 + x_move;
    x_marker_2 = 0.90 + x_move;
    x_text = 1.08 + x_move;
    y_rows = [-0.75, -1.35, -1.95] + 3.00 + y_move;
    swatch_x = 0.64 + x_move;
    swatch_width = 0.36;
    swatch_height = 0.36;

    patch(ax, [swatch_x, swatch_x + swatch_width, ...
        swatch_x + swatch_width, swatch_x], ...
        y_rows(1) + [-1, -1, 1, 1] * swatch_height / 2, 'r', ...
        'FaceAlpha', 0.14, 'EdgeColor', 'none');
    patch(ax, [swatch_x, swatch_x + swatch_width, ...
        swatch_x + swatch_width, swatch_x], ...
        y_rows(2) + [-1, -1, 1, 1] * swatch_height / 2, 'b', ...
        'FaceAlpha', 0.14, 'EdgeColor', 'none');
    patch(ax, [swatch_x, swatch_x + swatch_width, ...
        swatch_x + swatch_width, swatch_x], ...
        y_rows(3) + [-1, -1, 1, 1] * swatch_height / 2, [0.25, 0.25, 0.25], ...
        'FaceAlpha', 0.14, 'EdgeColor', 'none');

    plot(ax, x_marker_1, y_rows(1), marker_cfg.psi20.marker, 'Color', 'r', ...
        'MarkerFaceColor', 'r', ...
        'MarkerSize', marker_cfg.psi20.size, 'LineWidth', 1.2);
    plot(ax, x_marker_2, y_rows(1), marker_cfg.psi10.marker, 'Color', 'r', ...
        'MarkerFaceColor', 'none', ...
        'MarkerSize', marker_cfg.psi10.size, 'LineWidth', 1.2);
    text(ax, x_text, y_rows(1), '$R=20$ mm', ...
        'Interpreter', 'latex', 'FontName', 'Times New Roman', ...
        'FontSize', fs - 5, 'VerticalAlignment', 'middle');

    plot(ax, x_marker_1, y_rows(2), marker_cfg.psi20.marker, 'Color', 'b', ...
        'MarkerFaceColor', 'b', ...
        'MarkerSize', marker_cfg.psi20.size, 'LineWidth', 1.2);
    plot(ax, x_marker_2, y_rows(2), marker_cfg.psi10.marker, 'Color', 'b', ...
        'MarkerFaceColor', 'none', ...
        'MarkerSize', marker_cfg.psi10.size, 'LineWidth', 1.2);
    text(ax, x_text, y_rows(2), '$R=40$ mm', ...
        'Interpreter', 'latex', 'FontName', 'Times New Roman', ...
        'FontSize', fs - 5, 'VerticalAlignment', 'middle');

    plot(ax, (x_marker_1 + x_marker_2) / 2, y_rows(3), ...
        marker_cfg.reference.marker, 'Color', 'k', ...
        'MarkerFaceColor', 'k', ...
        'MarkerSize', marker_cfg.reference.size, 'LineWidth', 1.2);
    text(ax, x_text, y_rows(3), 'No channel', ...
        'FontName', 'Times New Roman', 'FontSize', fs - 5, ...
        'VerticalAlignment', 'middle');
end

function plot_pooled_errorbars(ax, stats, group_indices, color)
    lambda_1 = [];
    theta_1 = [];
    lambda_2 = [];
    theta_2 = [];
    for group_idx = group_indices
        lambda_1 = [lambda_1, stats.(sprintf('lam1_%d', group_idx))]; %#ok<AGROW>
        theta_1 = [theta_1, stats.(sprintf('deg1_%d', group_idx))]; %#ok<AGROW>
        lambda_2 = [lambda_2, stats.(sprintf('lam2_%d', group_idx))]; %#ok<AGROW>
        theta_2 = [theta_2, stats.(sprintf('deg2_%d', group_idx))]; %#ok<AGROW>
    end
    errorbar(ax, mean(lambda_1), mean(theta_1), std(theta_1), std(theta_1), ...
        std(lambda_1), std(lambda_1), 'Color', color, 'LineStyle', 'none', ...
        'Marker', 'none', 'LineWidth', 1.4, 'CapSize', 8);
    errorbar(ax, mean(lambda_2), mean(theta_2), std(theta_2), std(theta_2), ...
        std(lambda_2), std(lambda_2), 'Color', color, 'LineStyle', 'none', ...
        'Marker', 'none', 'LineWidth', 1.4, 'CapSize', 8);
end

function plot_centroid_connection(ax, stats, group_indices, color)
    lambda_centroid = zeros(1, 2);
    theta_centroid = zeros(1, 2);
    for spectral_group = 1:2
        lambda = [];
        theta = [];
        for group_idx = group_indices
            lambda = [lambda, ...
                stats.(sprintf('lam%d_%d', spectral_group, group_idx))]; %#ok<AGROW>
            theta = [theta, ...
                stats.(sprintf('deg%d_%d', spectral_group, group_idx))]; %#ok<AGROW>
        end
        valid = isfinite(lambda) & isfinite(theta);
        lambda_centroid(spectral_group) = mean(lambda(valid));
        theta_centroid(spectral_group) = mean(theta(valid));
    end
    plot(ax, lambda_centroid, theta_centroid, '--', 'Color', color, ...
        'LineWidth', 1.3);
end

function plot_group_shading(ax, stats, group_indices, color, cfg)
    for spectral_group = 1:2
        lambda = [];
        theta = [];
        for group_idx = group_indices
            lambda = [lambda, ...
                stats.(sprintf('lam%d_%d', spectral_group, group_idx))]; %#ok<AGROW>
            theta = [theta, ...
                stats.(sprintf('deg%d_%d', spectral_group, group_idx))]; %#ok<AGROW>
        end
        plot_cluster_envelope(ax, lambda, theta, color, cfg);
    end
end

function plot_cluster_envelope(ax, lambda, theta, color, cfg)
    valid = isfinite(lambda) & isfinite(theta);
    points = [lambda(valid)', theta(valid)'];
    points = unique(points, 'rows', 'stable');
    if size(points, 1) < 3
        return;
    end

    hull_idx = convhull(points(:, 1), points(:, 2));
    hull = points(hull_idx, :);
    center = mean(points, 1);
    hull = center + cfg.expansion * (hull - center);
    hull = smooth_closed_polygon(hull, cfg.smoothing_iterations);

    patch(ax, hull(:, 1), hull(:, 2), color, ...
        'FaceAlpha', cfg.alpha, 'EdgeColor', 'none');
end

function smoothed = smooth_closed_polygon(vertices, iterations)
    if isequal(vertices(1, :), vertices(end, :))
        vertices(end, :) = [];
    end
    smoothed = vertices;
    for iteration = 1:iterations
        following = circshift(smoothed, -1, 1);
        first_cut = 0.75 * smoothed + 0.25 * following;
        second_cut = 0.25 * smoothed + 0.75 * following;
        refined = zeros(2 * size(smoothed, 1), 2);
        refined(1:2:end, :) = first_cut;
        refined(2:2:end, :) = second_cut;
        smoothed = refined;
    end
    smoothed(end + 1, :) = smoothed(1, :);
end
