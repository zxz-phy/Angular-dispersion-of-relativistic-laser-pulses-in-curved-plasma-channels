clc;
clear;
close all;

%% Runtime options
if_export = true;  % true: save PNG/EPS/PDF; false: display only
axbc_left=0.1;
%% Geometrical-angle correction
% The raw angular calibration uses the discharge-off beam as its zero. Remove
% the wavelength-independent turning angle theta_0 = S/R so that theta is
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
if_center_NIRMIR = true;
if_cross = true;
center_cfg.split_wavelength_um = 2.5;
center_cfg.peak_fraction = 0.80;
center_cfg.line_style = '--';       % '-' for solid, '--' for dashed
center_cfg.size = 0.05;             % cross half-length relative to each axis span
center_cfg.line_width = 2;
center_cfg.color = 'k';
center_cfg.if_cross = if_cross;

%% Paths
script_dir = fileparts(mfilename('fullpath'));
figs_root = fileparts(fileparts(script_dir));
data_dir = fullfile(figs_root, 'data', 'Fig2');
output_dir = fullfile(figs_root, 'Figs', 'Fig2');
if if_export && ~isfolder(output_dir)
    mkdir(output_dir);
end

%% Figure style
fs = 20;
fg = figure(2);
set(fg, 'Color', 'w', 'Units', 'pixels', ...
    'Position', [200, 50, 1050, 830]);
% Explicit axes positions reserve a dedicated slot for the panel-(b) colorbar.
% This avoids tiledlayout's restriction on independently shifting the right column.
left_x = 0.095;
right_x = 0.515;
bottom_y = 0.105;
top_y = 0.565;
panel_width = 0.300;
panel_height = 0.335;

base_map = jet(1000);
n_fade = 40;
fade_to_white = [linspace(base_map(n_fade, 1), 1, n_fade)', ...
    linspace(base_map(n_fade, 2), 1, n_fade)', ...
    linspace(base_map(n_fade, 3), 1, n_fade)'];
intensity_map = [flipud(fade_to_white); base_map];

calibration = load(fullfile(data_dir, 'wavelength_calibration.mat'));
wavelength_table = calibration.wavelength_corrspond_concave_mirror(2:2:end, :);

%% (a) S = 2 mm, R = 20 mm
ax_a = axes('Parent', fg, 'Units', 'normalized', ...
    'Position', [left_x, top_y, panel_width, panel_height]);
image_a = read_pyro_spectrum(fullfile(data_dir, 'spectrum_R20mm.bgData'), ...
    4, 1.5, 20, 500, 350);
% 9.6354%
plot_spectrum(ax_a, image_a, 9.01,13.25, wavelength_table, ...
    [4, 23], theta_offset_a_deg, intensity_map, fs);
if if_center_NIRMIR
    plot_NIRMIR_centers(ax_a, image_a, 9.01, wavelength_table, center_cfg);
end
title(ax_a, '$S=2~\mathrm{mm},\ R=20~\mathrm{mm}$', ...
    'Interpreter', 'latex');
add_panel_label(ax_a, '(a)', fs);

%% (b) S = 2 mm, R = 40 mm
ax_b = axes('Parent', fg, 'Units', 'normalized', ...
    'Position', [right_x, top_y, panel_width, panel_height]);
image_b = read_pyro_spectrum(fullfile(data_dir, 'spectrum_R40mm.bgData'), ...
    16, 2, 2.5, 5000, 0);
plot_spectrum(ax_b, image_b, 1.53, 16.2, wavelength_table, ...
    [12, 23], theta_offset_b_deg, intensity_map, fs);
if if_center_NIRMIR
    plot_NIRMIR_centers(ax_b, image_b, 1.53, wavelength_table, center_cfg);
end
title(ax_b, '$S=2~\mathrm{mm},\ R=40~\mathrm{mm}$', ...
    'Interpreter', 'latex');
add_panel_label(ax_b, '(b)', fs);
cb = colorbar(ax_b, 'eastoutside');
cb.Label.String = 'Intensity (a. u.)';
cb.Label.FontName = 'Times New Roman';
cb.Label.FontSize = fs;
cb.Label.Rotation = 90;
drawnow;
cb.Position = [0.840, top_y, 0.018, panel_height];

%% (c) Without discharge
ax_c = axes('Parent', fg, 'Units', 'normalized', ...
    'Position', [left_x, bottom_y, panel_width, panel_height]);
image_c = read_pyro_spectrum(fullfile(data_dir, ...
    'spectrum_without_discharge.bgData'), 6, 2, 2.5, 500, 0);
plot_spectrum(ax_c, image_c, 9.7, 16.2, wavelength_table, ...
    [12, 23], theta_offset_c_deg, intensity_map, fs);
if if_center_NIRMIR
    plot_NIRMIR_centers(ax_c, image_c, 9.7, wavelength_table, center_cfg);
end
title(ax_c, 'Without discharge');
add_panel_label(ax_c, '(c)', fs);

%% (d) Statistical angular spectra
ax_d = axes('Parent', fg, 'Units', 'normalized', ...
    'Position', [right_x, bottom_y, panel_width, panel_height]);
stats = load(fullfile(data_dir, 'angular_dispersion_statistics.mat'));
stats = subtract_channel_turning_angle(stats, [1, 2], theta0_R20_deg);
stats = subtract_channel_turning_angle(stats, [3, 4], theta0_R40_deg);
hold(ax_d, 'on');

% Marker and uncertainty settings shared by the data and legend.
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

xlim(ax_d, [0.5, 5]);
ylim(ax_d, [-6.0, 0.5]);
% axd_values=[-6:2:0];
% yticklabels(ax_d, compose('%.1f', axd_values));

xlabel(ax_d, '$\lambda~(\mu\mathrm{m})$', 'Interpreter', 'latex');
ylabel(ax_d, '$\theta~(^{\circ})$', 'Interpreter', 'latex');
title(ax_d, 'Angular-dispersion statistics', 'FontSize', fs - 4);
set(ax_d, 'FontName', 'Times New Roman', 'FontSize', fs, ...
    'Box', 'on', 'TickDir', 'out');
add_grouped_legend(ax_d, fs, marker_cfg);
add_panel_label(ax_d, '(d)', fs);

%% Export
if if_export
    output_png = fullfile(output_dir, 'Fig2.png');
    output_eps = fullfile(output_dir, 'Fig2.eps');
    output_pdf = fullfile(output_dir, 'Fig2.pdf');
    exportgraphics(fg, output_png, 'Resolution', 300);
    print(fg, output_eps, '-depsc', '-painters');
    exportgraphics(fg, output_pdf, 'ContentType', 'vector');

    fprintf('Saved Fig. 2 to:\n%s\n%s\n%s\n', ...
        output_png, output_eps, output_pdf);
else
    fprintf('Fig. 2 displayed without export (if_export = false).\n');
end

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

function plot_spectrum(ax, spectrum, position_800nm, angle_origin, ...
        wavelength_table, y_limits, theta_offset_deg, map, fs)
    imagesc(ax, [0, 25.6], [0, 25.6], spectrum ./ max(spectrum(:)));
    wavelength_ticks = position_800nm + wavelength_table(:, 2);
    xticks(ax, wavelength_ticks);
    xticklabels(ax, string(wavelength_table(:, 1)));
    angle_values = (-2:0.5:2)+round(theta_offset_deg);
    angle_ticks = angle_origin -381 * tan(deg2rad(theta_offset_deg))+ 381 * tan(deg2rad(angle_values));
    
    if ax.Position(1)<0.13 && ax.Position(2)>0.55 %ax_a
        yticks(ax, angle_ticks(1:2:end));
        yticklabels(ax, compose('%.1f', angle_values(1:2:end)));
    else
        yticks(ax, angle_ticks);
        yticklabels(ax, compose('%.1f', angle_values));
    end
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

function plot_NIRMIR_centers(ax, spectrum, position_800nm, ...
        wavelength_table, cfg)
    % Split the detector image at the selected wavelength, retain the
    % brightest part of each band, and mark its intensity-weighted centroid.
    x_values = linspace(0, 25.6, size(spectrum, 2));
    y_values = linspace(0, 25.6, size(spectrum, 1));
    wavelength_positions = position_800nm + wavelength_table(:, 2);
    split_x = interp1(wavelength_table(:, 1), wavelength_positions, ...
        cfg.split_wavelength_um, 'linear', 'extrap');

    visible_x = xlim(ax);
    band_masks = {x_values >= visible_x(1) & x_values < split_x, ...
        x_values >= split_x & x_values <= visible_x(2)};
    hold_state = ishold(ax);
    hold(ax, 'on');
    for band_idx = 1:numel(band_masks)
        column_mask = band_masks{band_idx};
        band = spectrum(:, column_mask);
        if isempty(band) || ~any(isfinite(band(:)))
            continue;
        end

        band_max = max(band(:), [], 'omitnan');
        if ~(isfinite(band_max) && band_max > 0)
            continue;
        end

        retained = isfinite(band) & band >= cfg.peak_fraction * band_max;
        weights = band;
        weights(~retained) = 0;
        weight_sum = sum(weights(:));
        if weight_sum <= 0
            continue;
        end

        [x_grid, y_grid] = meshgrid(x_values(column_mask), y_values);
        center_x = sum(x_grid(:) .* weights(:)) / weight_sum;
        center_y = sum(y_grid(:) .* weights(:)) / weight_sum;
        draw_center_marker(ax, center_x, center_y, cfg);
    end
    if ~hold_state
        hold(ax, 'off');
    end
end

function draw_center_marker(ax, center_x, center_y, cfg)
    % Select either a centered cross or projections to the coordinate axes.
    x_limits = xlim(ax);
    y_limits = ylim(ax);
    if cfg.if_cross
        half_width = cfg.size * diff(x_limits);
        half_height = cfg.size * diff(y_limits);
        plot_symmetric_cross(ax, center_x, center_y, ...
            half_width, half_height, cfg);
    else
        horizontal_x = [x_limits(1), center_x];
        horizontal_y = [center_y, center_y];
        vertical_x = [center_x, center_x];
        vertical_y = [y_limits(1), center_y];
        plot_centered_half_line(ax, horizontal_x, horizontal_y, cfg);
        plot_centered_half_line(ax, vertical_x, vertical_y, cfg);
    end
end

function plot_symmetric_cross(ax, center_x, center_y, ...
        half_width, half_height, cfg)
    if strcmp(cfg.line_style, '--')
        % MATLAB anchors its dash pattern at a line endpoint. Explicitly
        % mirrored segments keep a short dashed cross symmetric about the
        % centroid and ensure that the dashes remain visually distinct.
        segments = [-1.00, -0.48; -0.18, 0.18; 0.48, 1.00];
        for segment_idx = 1:size(segments, 1)
            x_segment = center_x + half_width * segments(segment_idx, :);
            y_segment = center_y + half_height * segments(segment_idx, :);
            plot_marker_line(ax, x_segment, [center_y, center_y], '-', cfg);
            plot_marker_line(ax, [center_x, center_x], y_segment, '-', cfg);
        end
    else
        plot_marker_line(ax, center_x + [-half_width, half_width], ...
            [center_y, center_y], cfg.line_style, cfg);
        plot_marker_line(ax, [center_x, center_x], ...
            center_y + [-half_height, half_height], cfg.line_style, cfg);
    end
end

function plot_centered_half_line(ax, x_values, y_values, cfg)
    plot_marker_line(ax, x_values, y_values, cfg.line_style, cfg);
end

function plot_marker_line(ax, x_values, y_values, line_style, cfg)
    plot(ax, x_values, y_values, 'LineStyle', line_style, 'Color', cfg.color, ...
        'LineWidth', cfg.line_width, 'HandleVisibility', 'off');
end

function stats = subtract_channel_turning_angle(stats, group_indices, theta0_deg)
    % Shift both spectral clusters for every pressure condition belonging
    % to the same channel geometry. The discharge-off reference is retained
    % in its original reference coordinate.
    for group_idx = group_indices
        for spectral_group = 1:2
            field_name = sprintf('deg%d_%d', spectral_group, group_idx);
            stats.(field_name) = stats.(field_name) - theta0_deg;
        end
    end
end

function add_panel_label(ax, label, fs)
    text(ax, -0.160, 1.035, label, 'Units', 'normalized', ...
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
    % Keep the legend inside the data limits so its frame is not clipped.
    x_left = 3.00;
    y_bottom = -2.40;
    legend_width = 1.75;
    legend_height = 1.55;
    rectangle(ax, 'Position', [x_left, y_bottom, legend_width, legend_height], ...
        'FaceColor', 'w', 'EdgeColor', 'k', 'LineWidth', 0.8);
    x_marker_1 = x_left + 0.15;
    x_marker_2 = x_left + 0.33;
    x_text = x_left + 0.50;
    y_rows = y_bottom + [1.22, 0.77, 0.32];
    swatch_x = x_left + 0.07;
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
        'FontSize', fs - 7, 'VerticalAlignment', 'middle');

    plot(ax, x_marker_1, y_rows(2), marker_cfg.psi20.marker, 'Color', 'b', ...
        'MarkerFaceColor', 'b', ...
        'MarkerSize', marker_cfg.psi20.size, 'LineWidth', 1.2);
    plot(ax, x_marker_2, y_rows(2), marker_cfg.psi10.marker, 'Color', 'b', ...
        'MarkerFaceColor', 'none', ...
        'MarkerSize', marker_cfg.psi10.size, 'LineWidth', 1.2);
    text(ax, x_text, y_rows(2), '$R=40$ mm', ...
        'Interpreter', 'latex', 'FontName', 'Times New Roman', ...
        'FontSize', fs - 7, 'VerticalAlignment', 'middle');

    plot(ax, (x_marker_1 + x_marker_2) / 2, y_rows(3), ...
        marker_cfg.reference.marker, 'Color', 'k', ...
        'MarkerFaceColor', 'k', ...
        'MarkerSize', marker_cfg.reference.size, 'LineWidth', 1.2);
    text(ax, x_text, y_rows(3), 'No channel', ...
        'FontName', 'Times New Roman', 'FontSize', fs - 7, ...
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
