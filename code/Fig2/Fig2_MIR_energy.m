clc;
clear;
close all;

%% Paths
script_dir = fileparts(mfilename('fullpath'));
figs_root = fileparts(fileparts(script_dir));
data_dir = fullfile(figs_root, 'data', 'Fig2');

% measured relative to the corresponding curved-channel axis.
S_mm = 2;
R20_mm = 20;
R40_mm = 40;
theta0_R20_deg = rad2deg(S_mm / R20_mm);
% theta0_R20_deg = 0;
theta0_R40_deg = rad2deg(S_mm / R40_mm);
% theta0_R40_deg = 0;
theta_offset_a_deg = -theta0_R20_deg;
theta_offset_b_deg = -theta0_R40_deg;
theta_offset_c_deg = -theta0_R40_deg;

%% Spectral-analysis settings
nir_range_um = [0.8, 2.0];
mir_range_um = [2.5, 5.0];
mir_energy_range_um = [3.0, 4.5];

%% Figure style
fs = 20;
fg = figure(2);
set(fg, 'Color', 'w', 'Units', 'pixels', ...
    'Position', [200, 50, 1050, 830]);
layout = tiledlayout(fg, 2, 2, 'TileSpacing', 'compact', ...
    'Padding', 'compact');

base_map = jet(1000);
n_fade = 40;
fade_to_white = [linspace(base_map(n_fade, 1), 1, n_fade)', ...
    linspace(base_map(n_fade, 2), 1, n_fade)', ...
    linspace(base_map(n_fade, 3), 1, n_fade)'];
intensity_map = [flipud(fade_to_white); base_map];

calibration = load(fullfile(data_dir, 'wavelength_calibration.mat'));
wavelength_calibration = calibration.wavelength_corrspond_concave_mirror;
wavelength_table = calibration.wavelength_corrspond_concave_mirror(2:2:end, :);

%% (a) S = 2 mm, R = 20 mm
ax_a = nexttile(layout, 1);
image_a = read_pyro_spectrum(fullfile(data_dir, 'spectrum_R20mm.bgData'), ...
    4, 1.5, 20, 500, 350);
%9.6354
plot_spectrum(ax_a, image_a, 9.01, 13.25, wavelength_table, ...
    [4, 23], intensity_map, fs);
title(ax_a, '$S=2~\mathrm{mm},\ R=20~\mathrm{mm}$', ...
    'Interpreter', 'latex');
add_panel_label(ax_a, '(a)', fs);
result_a = annotate_spectral_centroids(ax_a, image_a, 9.01, ...
    wavelength_calibration, nir_range_um, mir_range_um, true, fs);
mir_signal_a = estimate_band_signal(fullfile(data_dir, 'spectrum_R20mm.bgData'), ...
    4, 1.5, 9.01, wavelength_calibration, mir_energy_range_um);
nir_signal_a = estimate_band_signal(fullfile(data_dir, 'spectrum_R20mm.bgData'), ...
    4, 1.5, 9.01, wavelength_calibration, nir_range_um);

%% (b) S = 2 mm, R = 40 mm
ax_b = nexttile(layout, 2);
image_b = read_pyro_spectrum(fullfile(data_dir, 'spectrum_R40mm.bgData'), ...
    16, 2, 2.5, 5000, 0);
plot_spectrum(ax_b, image_b, 1.53, 16.2, wavelength_table, ...
    [12, 23], intensity_map, fs);
title(ax_b, '$S=2~\mathrm{mm},\ R=40~\mathrm{mm}$', ...
    'Interpreter', 'latex');
add_panel_label(ax_b, '(b)', fs);
result_b = annotate_spectral_centroids(ax_b, image_b, 1.53, ...
    wavelength_calibration, nir_range_um, mir_range_um, true, fs);
mir_signal_b = estimate_band_signal(fullfile(data_dir, 'spectrum_R40mm.bgData'), ...
    16, 2, 1.53, wavelength_calibration, mir_energy_range_um);
nir_signal_b = estimate_band_signal(fullfile(data_dir, 'spectrum_R40mm.bgData'), ...
    16, 2, 1.53, wavelength_calibration, nir_range_um);
cb = colorbar(ax_b);
cb.Label.String = 'Intensity (a. u.)';
cb.Label.FontName = 'Times New Roman';
cb.Label.FontSize = fs;

%% (c) Without discharge
ax_c = nexttile(layout, 3);
image_c = read_pyro_spectrum(fullfile(data_dir, ...
    'spectrum_without_discharge.bgData'), 6, 2, 2.5, 500, 0);
plot_spectrum(ax_c, image_c, 9.7, 16.2, wavelength_table, ...
    [12, 23], intensity_map, fs);
title(ax_c, 'Without discharge');
add_panel_label(ax_c, '(c)', fs);
result_c = annotate_spectral_centroids(ax_c, image_c, 9.7, ...
    wavelength_calibration, nir_range_um, mir_range_um, false, fs);
mir_signal_c = estimate_band_signal(fullfile(data_dir, ...
    'spectrum_without_discharge.bgData'), 6, 2, 9.7, wavelength_calibration, mir_energy_range_um);
nir_signal_c = estimate_band_signal(fullfile(data_dir, ...
    'spectrum_without_discharge.bgData'), 6, 2, 9.7, wavelength_calibration, nir_range_um);

%% (d) Statistical angular spectra
ax_d = nexttile(layout, 4);
stats = load(fullfile(data_dir, 'angular_dispersion_statistics.mat'));
hold(ax_d, 'on');

% Marker and uncertainty settings shared by the data and legend.
if_errorbar = false;
if_shading = true;
if_connection = true;
shading_cfg.alpha = 0.2;
shading_cfg.expansion = 1.08;
shading_cfg.smoothing_iterations = 5;
psi10_size = 5;
psi20_size = 5;
reference_size = 5;
marker_cfg.psi20.marker = 'o';
marker_cfg.psi20.size = psi20_size;
marker_cfg.psi20.filled = true;
marker_cfg.psi10.marker = 'o';
marker_cfg.psi10.size = psi10_size;
marker_cfg.psi10.filled = false;
marker_cfg.reference.marker = 'o';
marker_cfg.reference.size = reference_size;
marker_cfg.reference.filled = true;

if if_shading
    % Shade the NIR and MIR clusters separately for each condition.
    plot_group_shading(ax_d, stats, [3, 4], [0.20, 0.35, 1.00], shading_cfg);
    plot_group_shading(ax_d, stats, [1, 2], [1.00, 0.20, 0.20], shading_cfg);
    plot_group_shading(ax_d, stats, 5, [0.25, 0.25, 0.25], shading_cfg);
end

if if_connection
    % Connect the centroids of the short- and long-wavelength clusters.
    plot_centroid_connection(ax_d, stats, [3, 4], 'b');
    plot_centroid_connection(ax_d, stats, [1, 2], 'r');
    plot_centroid_connection(ax_d, stats, 5, 'k');
end

if if_errorbar
    % Draw uncertainties below the raw markers. One pooled uncertainty
    % series is used for each channel condition, irrespective of pressure.
    plot_pooled_errorbars(ax_d, stats, [3, 4], 'b');
    plot_pooled_errorbars(ax_d, stats, [1, 2], 'r');
    plot_pooled_errorbars(ax_d, stats, 5, 'k');
end

% Draw black and blue first so the less numerous red points remain visible.
plot_pair(ax_d, stats.lam1_5, stats.deg1_5, ...
    stats.lam2_5, stats.deg2_5, marker_cfg.reference, 'k');
plot_pair(ax_d, stats.lam1_3, stats.deg1_3, ...
    stats.lam2_3, stats.deg2_3, marker_cfg.psi20, 'b');
plot_pair(ax_d, stats.lam1_4, stats.deg1_4, ...
    stats.lam2_4, stats.deg2_4, marker_cfg.psi10, 'b');
plot_pair(ax_d, stats.lam1_1, stats.deg1_1, ...
    stats.lam2_1, stats.deg2_1, marker_cfg.psi20, 'r');
plot_pair(ax_d, stats.lam1_2, stats.deg1_2, ...
    stats.lam2_2, stats.deg2_2, marker_cfg.psi10, 'r');

% Wavelength statistics displayed in panel (d). Groups 1--2 and 3--4 pool
% the two pressure cases for the same channel geometry.
wavelength_stats.r20 = summarize_case_wavelengths(stats, [1, 2]);
wavelength_stats.r40 = summarize_case_wavelengths(stats, [3, 4]);
wavelength_stats.off = summarize_case_wavelengths(stats, 5);

angle_stats.r20 = summarize_case_angles(stats, [1, 2],theta_offset_a_deg);
angle_stats.r40 = summarize_case_angles(stats, [3, 4],theta_offset_b_deg);
angle_stats.off = summarize_case_angles(stats, 5,theta_offset_c_deg);

xlim(ax_d, [0.5, 5]);
ylim(ax_d, [-1.2, 1.0]);
xlabel(ax_d, '$\lambda~(\mu\mathrm{m})$', 'Interpreter', 'latex');
ylabel(ax_d, '$\theta~(^{\circ})$', 'Interpreter', 'latex');
title(ax_d, 'Angular spectra');
set(ax_d, 'FontName', 'Times New Roman', 'FontSize', fs, ...
    'Box', 'on', 'TickDir', 'out');
add_grouped_legend(ax_d, fs, marker_cfg);
display_wavelength_statistics(ax_d, wavelength_stats, fs);
add_panel_label(ax_d, '(d)', fs);

%% Report numerical results (this analysis script intentionally does not export)
print_spectral_result('(a) S=2 mm, R=20 mm', result_a);
print_spectral_result('(b) S=2 mm, R=40 mm', result_b);
print_spectral_result('(c) Without discharge', result_c);
print_signal_comparison(mir_signal_a, mir_signal_b, mir_signal_c, ...
    nir_signal_a, nir_signal_b, nir_signal_c);
print_case_wavelength_statistics(wavelength_stats);
print_case_angle_statistics(angle_stats)

function spectrum = read_pyro_spectrum(file_path, shot_index, ...
        rotation_deg, threshold, min_area, zero_columns)
    dataset = sprintf('/BG_DATA/%d/DATA', shot_index);
    height_path = sprintf('/BG_DATA/%d/RAWFRAME/HEIGHT', shot_index);
    width_path = sprintf('/BG_DATA/%d/RAWFRAME/WIDTH', shot_index);
    raw = h5read(file_path, dataset);
    height = h5read(file_path, height_path);
    width = h5read(file_path, width_path);
    image = reshape(raw, width, height);
    image = fliplr(double(image)' ./ 2^19);
    image = interp2(image, 3);
    image = imrotate(image, rotation_deg);
    if rotation_deg == 1.5
        gaussian_sigma = 3;
    else
        gaussian_sigma = 2;
    end
    image = imgaussfilt(image, gaussian_sigma);
    mask = imbinarize(image, threshold);
    mask = bwareaopen(mask, min_area);
    if zero_columns > 0
        mask(:, 1:min(zero_columns, size(mask, 2))) = 0;
    end
    spectrum = image .* mask;
end

function signal = estimate_band_signal(file_path, shot_index, rotation_deg, ...
        position_800nm, wavelength_calibration, mir_range_um)
    % Use the unnormalized camera frame for relative signal integration.
    image = read_pyro_frame(file_path, shot_index, rotation_deg);
    x_axis = linspace(0, 25.6, size(image, 2));
    calibration_x = position_800nm + wavelength_calibration(:, 2);
    position_coeff = polyfit(wavelength_calibration(:, 1), calibration_x, 2);
    lambda_axis = invert_quadratic_position(position_coeff, x_axis, ...
        [min(wavelength_calibration(:, 1)), max(wavelength_calibration(:, 1))]);
    region = image(:, lambda_axis >= mir_range_um(1) & lambda_axis < mir_range_um(2));

    % Estimate the local camera background from the lowest-intensity pixels so
    % broad NIR features are not absorbed into the background estimate.
    background_samples = region(region <= prctile(region(:), 20));
    signal.background = median(background_samples);
    signal.noise = 1.4826 * mad(background_samples, 1);
    mask = region > signal.background + 5 * signal.noise;
    mask = bwareaopen(mask, 50);
    signal.counts = sum(region(mask) - signal.background);
    signal.n_pixels = nnz(mask);
end

function image = read_pyro_frame(file_path, shot_index, rotation_deg)
    dataset = sprintf('/BG_DATA/%d/DATA', shot_index);
    height_path = sprintf('/BG_DATA/%d/RAWFRAME/HEIGHT', shot_index);
    width_path = sprintf('/BG_DATA/%d/RAWFRAME/WIDTH', shot_index);
    raw = h5read(file_path, dataset);
    height = h5read(file_path, height_path);
    width = h5read(file_path, width_path);
    image = reshape(raw, width, height);
    image = fliplr(double(image)' ./ 2^19);
    image = interp2(image, 3);
    image = imrotate(image, rotation_deg);
    if rotation_deg == 1.5
        gaussian_sigma = 3;
    else
        gaussian_sigma = 2;
    end
    image = imgaussfilt(image, gaussian_sigma);
end

function plot_spectrum(ax, spectrum, position_800nm, angle_origin, ...
        wavelength_table, y_limits, map, fs)
    imagesc(ax, [0, 25.6], [0, 25.6], spectrum ./ max(spectrum(:)));
    wavelength_ticks = position_800nm + wavelength_table(:, 2);
    xticks(ax, wavelength_ticks);
    xticklabels(ax, string(wavelength_table(:, 1)));
    angle_values = -2:0.5:2;
    angle_ticks = angle_origin + 381 * tan(deg2rad(angle_values));
    yticks(ax, angle_ticks);
    yticklabels(ax, string(angle_values));
    xlim(ax, [position_800nm - 0.1, wavelength_ticks(5)]);
    ylim(ax, y_limits);
    caxis(ax, [0, 1]);
    colormap(ax, map);
    set(ax, 'YDir', 'normal', 'FontName', 'Times New Roman', ...
        'FontSize', fs, 'Box', 'on', 'TickDir', 'out', ...
        'TickLength', [0.02, 0.025]);
    xlabel(ax, '$\lambda~(\mu\mathrm{m})$', 'Interpreter', 'latex');
    ylabel(ax, '$\theta~(^{\circ})$', 'Interpreter', 'latex');
end

function result = annotate_spectral_centroids(ax, spectrum, position_800nm, ...
        wavelength_calibration, nir_range_um, mir_range_um, ...
        show_mir, fs)
    % Approximate detector position as a quadratic function of wavelength,
    % x(lambda), and invert it to obtain a continuous wavelength axis.
    calibration_x = position_800nm + wavelength_calibration(:, 2);
    calibration_lambda = wavelength_calibration(:, 1);
    position_coeff = polyfit(calibration_lambda, calibration_x, 2);
    calibration_range = [min(calibration_lambda), max(calibration_lambda)];

    n_rows = size(spectrum, 1);
    n_cols = size(spectrum, 2);
    x_axis = linspace(0, 25.6, n_cols);
    y_axis = linspace(0, 25.6, n_rows);
    lambda_axis = invert_quadratic_position(position_coeff, x_axis, ...
        calibration_range);

    result.nir = spectral_region_centroid(spectrum, x_axis, y_axis, ...
        lambda_axis, nir_range_um, position_coeff, calibration_range);
    result.mir = spectral_region_centroid(spectrum, x_axis, y_axis, ...
        lambda_axis, mir_range_um, position_coeff, calibration_range);
    result.nir_to_mir_area_ratio = NaN;
    if show_mir && result.mir.integral > 0
        result.nir_to_mir_area_ratio = result.nir.integral / result.mir.integral;
    end

    hold(ax, 'on');
    nir_color = [0.00, 0.35, 0.90];
    mir_color = [0.85, 0.05, 0.15];
    energy_gray = [0.42, 0.42, 0.42];

    plot_centroid_marker(ax, result.nir, nir_color);
    text(ax, result.nir.x, result.nir.y, ...
        sprintf('  NIR: %.2f \\mum', result.nir.lambda_um), ...
        'FontName', 'Times New Roman', 'FontSize', fs - 5, ...
        'FontWeight', 'bold', 'Color', nir_color, ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', ...
        'Clipping', 'on');

    if show_mir && result.mir.integral > 0
        plot_centroid_marker(ax, result.mir, mir_color);
        text(ax, result.mir.x, result.mir.y, ...
            sprintf('  MIR: %.2f \\mum', result.mir.lambda_um), ...
            'FontName', 'Times New Roman', 'FontSize', fs - 5, ...
            'FontWeight', 'bold', 'Color', mir_color, ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', ...
            'Clipping', 'on');
        text(ax, 0.97, 0.05, ...
            sprintf('NIR/MIR area ratio: %.1f', result.nir_to_mir_area_ratio), ...
            'Units', 'normalized', 'FontName', 'Times New Roman', ...
            'FontSize', fs - 5, 'Color', energy_gray, ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
    end
end

function centroid = spectral_region_centroid(spectrum, x_axis, y_axis, ...
        lambda_axis, lambda_range, position_coeff, calibration_range)
    column_mask = lambda_axis >= lambda_range(1) & ...
        lambda_axis < lambda_range(2);
    region = spectrum(:, column_mask);
    region(region<0.9*max(region(:))) = 0;
    total = sum(region(:));
    centroid.integral = total;
    centroid.x = NaN;
    centroid.y = NaN;
    centroid.lambda_um = NaN;
    if total <= 0 || ~any(column_mask)
        return;
    end

    x_region = x_axis(column_mask);
    column_signal = sum(region, 1);
    row_signal = sum(region, 2);
    centroid.x = sum(x_region .* column_signal) / total;
    centroid.y = sum(y_axis(:) .* row_signal) / total;
    centroid.lambda_um = invert_quadratic_position(position_coeff, ...
        centroid.x, calibration_range);
end

function lambda = invert_quadratic_position(position_coeff, x, valid_range)
    a = position_coeff(1);
    b = position_coeff(2);
    c = position_coeff(3) - x;
    discriminant = max(b.^2 - 4 .* a .* c, 0);
    root_1 = (-b + sqrt(discriminant)) ./ (2 .* a);
    root_2 = (-b - sqrt(discriminant)) ./ (2 .* a);
    lambda = root_1;
    use_root_2 = (root_1 < valid_range(1) | root_1 > valid_range(2)) & ...
        root_2 >= valid_range(1) & root_2 <= valid_range(2);
    lambda(use_root_2) = root_2(use_root_2);
end

function plot_centroid_marker(ax, centroid, color)
    if ~isfinite(centroid.x) || ~isfinite(centroid.y)
        return;
    end
    plot(ax, centroid.x, centroid.y, 'o', 'MarkerSize', 9, ...
        'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
    plot(ax, centroid.x, centroid.y, '+', 'MarkerSize', 8, ...
        'Color', color, 'LineWidth', 2);
end

function print_spectral_result(label, result)
    fprintf('\n%s\n', label);
    fprintf('  NIR centroid wavelength (quadratic calibration): %.4f um\n', ...
        result.nir.lambda_um);
    if isfinite(result.mir.lambda_um)
        fprintf('  MIR centroid wavelength (quadratic calibration): %.4f um\n', ...
            result.mir.lambda_um);
    end
    if isfinite(result.nir_to_mir_area_ratio)
        fprintf('  NIR/MIR area ratio from centroid masks: %.3f\n', ...
            result.nir_to_mir_area_ratio);
    else
        fprintf('  NIR/MIR area ratio: not evaluated (no reliable MIR reference).\n');
    end
end

function print_signal_comparison(mir_a, mir_b, mir_c, nir_a, nir_b, nir_c)
    fprintf('\nRelative 3.0--4.5 um MIR-band camera counts (background-subtracted)\n');
    fprintf('  R=20 mm: %.6g counts (%d pixels)\n', mir_a.counts, mir_a.n_pixels);
    fprintf('  R=40 mm: %.6g counts (%d pixels)\n', mir_b.counts, mir_b.n_pixels);
    fprintf('  Without discharge: %.6g counts (%d pixels)\n', mir_c.counts, mir_c.n_pixels);
    fprintf('  Without discharge / R=20 mm: %.4f\n', mir_c.counts / mir_a.counts);
    fprintf('  Without discharge / R=40 mm: %.4f\n', mir_c.counts / mir_b.counts);

    eta_a = mir_a.counts / nir_a.counts;
    eta_b = mir_b.counts / nir_b.counts;
    eta_c = mir_c.counts / nir_c.counts;
    fprintf('\nMIR/NIR band-count ratios\n');
    fprintf('  R=20 mm: %.6g\n', eta_a);
    fprintf('  R=40 mm: %.6g\n', eta_b);
    fprintf('  Without discharge: %.6g\n', eta_c);
    fprintf('  Without discharge / R=20 mm: %.4f\n', eta_c / eta_a);
    fprintf('  Without discharge / R=40 mm: %.4f\n', eta_c / eta_b);
end

function result = summarize_case_wavelengths(stats, group_indices)
    nir = [];
    mir = [];
    for group_idx = group_indices
        nir = [nir, stats.(sprintf('lam1_%d', group_idx))]; %#ok<AGROW>
        mir = [mir, stats.(sprintf('lam2_%d', group_idx))]; %#ok<AGROW>
    end
    nir = nir(isfinite(nir));
    mir = mir(isfinite(mir));
    result.nir_mean = mean(nir);
    result.nir_std = std(nir);
    result.mir_mean = mean(mir);
    result.mir_std = std(mir);
    result.nir_n = numel(nir);
    result.mir_n = numel(mir);
end

function result = summarize_case_angles(stats, group_indices,theta_offset_deg)
    nir_angle = [];
    mir_angle = [];
    for group_idx = group_indices
        nir_angle = [nir_angle, stats.(sprintf('deg1_%d', group_idx))+theta_offset_deg]; %#ok<AGROW>
        mir_angle = [mir_angle, stats.(sprintf('deg2_%d', group_idx))+theta_offset_deg]; %#ok<AGROW>
    end
    nir_angle = nir_angle(isfinite(nir_angle));
    mir_angle = mir_angle(isfinite(mir_angle));
    result.nir_angle_mean = mean(nir_angle);
    result.nir_angle_std = std(nir_angle);
    result.mir_angle_mean = mean(mir_angle);
    result.mir_angle_std = std(mir_angle);
    result.nir_angle_n = numel(nir_angle);
    result.mir_angle_n = numel(mir_angle);
end

function display_wavelength_statistics(ax, wavelength_stats, fs)
    entries = { ...
        wavelength_stats.r20, [0.85, 0.05, 0.05], '$R=20$ mm'; ...
        wavelength_stats.r40, [0.05, 0.20, 0.85], '$R=40$ mm'; ...
        wavelength_stats.off, [0.30, 0.30, 0.30], 'Discharge off'};
    y_positions = [0.94, 0.87, 0.80];
    for idx = 1:size(entries, 1)
        item = entries{idx, 1};
        label = entries{idx, 3};
        line = sprintf(['%s: $\\lambda_{\\rm NIR}=%.2f\\pm%.2f$, ', ...
            '$\\lambda_{\\rm MIR}=%.2f\\pm%.2f~\\mu$m'], ...
            label, item.nir_mean, item.nir_std, item.mir_mean, item.mir_std);
        text(ax, 0.98, y_positions(idx), line, 'Units', 'normalized', ...
            'Interpreter', 'latex', 'FontName', 'Times New Roman', ...
            'FontSize', fs - 8, 'Color', entries{idx, 2}, ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
            'BackgroundColor', 'w', 'Margin', 1.5);
    end
end

function print_case_angle_statistics(angle_states)
    fprintf('\nPanel (d) pooled angle statistics (mean +/- std, um)\n');
    print_one_case_angle('S=2 mm, R=20 mm', angle_states.r20);
    print_one_case_angle('S=2 mm, R=40 mm', angle_states.r40);
    print_one_case_angle('Without discharge', angle_states.off);
end

function print_case_wavelength_statistics(wavelength_stats)
    fprintf('\nPanel (d) pooled wavelength statistics (mean +/- std, um)\n');
    print_one_case('S=2 mm, R=20 mm', wavelength_stats.r20);
    print_one_case('S=2 mm, R=40 mm', wavelength_stats.r40);
    print_one_case('Without discharge', wavelength_stats.off);
end

function print_one_case(label, item)
    fprintf('  %s: NIR %.1f +/- %.1f (n=%d), MIR %.1f +/- %.1f (n=%d)\n', ...
        label, item.nir_mean, item.nir_std, item.nir_n, ...
        item.mir_mean, item.mir_std, item.mir_n);
end

function print_one_case_angle(label, item)
    fprintf('  %s: NIR %.1f +/- %.1f (n=%d), MIR %.1f +/- %.1f (n=%d)\n', ...
        label, item.nir_angle_mean, item.nir_angle_std, item.nir_angle_n, ...
        item.mir_angle_mean, item.mir_angle_std, item.mir_angle_n);
end

function add_panel_label(ax, label, fs)
    text(ax, -0.095, 1.035, label, 'Units', 'normalized', ...
        'FontName', 'Times New Roman', 'FontSize', fs + 4, ...
        'FontWeight', 'bold', 'VerticalAlignment', 'bottom', ...
        'HorizontalAlignment', 'left', 'Color', 'k', 'Clipping', 'off');
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
    % Three-row boxed legend in the lower-left corner.
    rectangle(ax, 'Position', [0.58, -1.145, 3, 0.58], ...
        'FaceColor', 'w', 'EdgeColor', 'k', 'LineWidth', 0.8);
    x_marker_1 = 0.72;
    x_marker_2 = 0.90;
    x_text = 1.08;
    y_rows = [-0.47, -0.64, -0.81]-0.22;
    swatch_x = 0.68;
    swatch_width = 0.36;
    swatch_height = 0.105;

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
    text(ax, x_text, y_rows(1), '$S=2$ mm, $R=20$ mm', ...
        'Interpreter', 'latex', 'FontName', 'Times New Roman', ...
        'FontSize', fs - 5, 'VerticalAlignment', 'middle');

    plot(ax, x_marker_1, y_rows(2), marker_cfg.psi20.marker, 'Color', 'b', ...
        'MarkerFaceColor', 'b', ...
        'MarkerSize', marker_cfg.psi20.size, 'LineWidth', 1.2);
    plot(ax, x_marker_2, y_rows(2), marker_cfg.psi10.marker, 'Color', 'b', ...
        'MarkerFaceColor', 'none', ...
        'MarkerSize', marker_cfg.psi10.size, 'LineWidth', 1.2);
    text(ax, x_text, y_rows(2), '$S=2$ mm, $R=40$ mm', ...
        'Interpreter', 'latex', 'FontName', 'Times New Roman', ...
        'FontSize', fs - 5, 'VerticalAlignment', 'middle');

    plot(ax, (x_marker_1 + x_marker_2) / 2, y_rows(3), ...
        marker_cfg.reference.marker, 'Color', 'k', ...
        'MarkerFaceColor', 'k', ...
        'MarkerSize', marker_cfg.reference.size, 'LineWidth', 1.2);
    text(ax, x_text, y_rows(3), 'Without discharge', ...
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
    % Keep the short- and long-wavelength clusters as separate envelopes.
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
    % Chaikin corner cutting rounds a closed polygon without spline overshoot.
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
