clc;
clear;
close all;

%% Paths and source files
script_dir = fileparts(mfilename('fullpath'));
figs_root = fileparts(fileparts(script_dir));
data_dir = fullfile(figs_root, 'data', 'Fig3');
output_dir = fullfile(figs_root, 'Figs', 'Fig3');
if ~isfolder(output_dir)
    mkdir(output_dir);
end

file_no_filter = fullfile(data_dir, ...
    '20250718_OAP+透镜pyro成像_18kV20psi.bgData');
file_no_filter_23psi = fullfile(data_dir, ...
    '20250718_OAP+透镜pyro成像_18kV23psi.bgData');
file_no_filter_34500fs = fullfile(data_dir, ...
    '20250718_OAP+透镜pyro成像_34500fs18kV20psi.bgData');
file_with_filter = fullfile(data_dir, ...
    '20250718_OAP+透镜pyro成像_外部加硅片_34500fs18kV20psi_真红外信号.bgData');
file_with_filter_2 = fullfile(data_dir, ...
    '20250718_OAP+透镜pyro成像_外部加硅片_34500fs18kV20psi_真红外信号_微弱偏下.bgData');
file_with_filter_3 = fullfile(data_dir, ...
    '20250718_OAP+透镜pyro成像_外部加硅片_34500fs18kV20psi_真红外信号_微弱偏下2.bgData');
file_no_discharge = fullfile(data_dir, ...
    '20250718_OAP+透镜pyro成像_34500fs不放电20psi.bgData');
file_no_discharge_10psi = fullfile(data_dir, ...
    '20250718_OAP+透镜pyro成像_不放电10psi.bgData');
file_no_discharge_with_filter_7j = fullfile(data_dir, ...
    '20250718_OAP+透镜pyro成像_外部加硅片_34500fs不放电20psi_7J_红外信号.bgData');
file_no_discharge_with_filter = fullfile(data_dir, ...
    '20250718_OAP+透镜pyro成像_外部加硅片_34500fs不放电20psi_红外信号.bgData');

%% Processing settings: no thresholding or connected-component masking
rotation_deg = 0;
gaussian_sigma = 1.0;
% Hard cuts are used only for centroid extraction. They do not alter the
% lightly processed images displayed in panels (a)-(c).
centroid_cut.no_filter = 0.45;
centroid_cut.with_filter = 0.40;
centroid_cut.no_discharge = 0.35;
centroid_cut.no_discharge_with_filter = 0.20;
centroid_min_area = 50;
centroid_min_area_weak = 15;
centroid_locator_sigma = 12;
centroid_roi_half_width = 0.18;
fs = 16;

no_filter{1} = read_bgdata_shot(file_no_filter, 1, rotation_deg, gaussian_sigma);
no_filter{2} = read_bgdata_shot(file_no_filter_23psi, 1, rotation_deg, gaussian_sigma);
no_filter{3} = read_bgdata_shot(file_no_filter_34500fs, 1, rotation_deg, gaussian_sigma);
with_filter{1} = read_bgdata_shot(file_with_filter, 1, rotation_deg, gaussian_sigma);
with_filter{2} = read_bgdata_shot(file_with_filter_2, 1, rotation_deg, gaussian_sigma);
with_filter{3} = read_bgdata_shot(file_with_filter_3, 1, rotation_deg, gaussian_sigma);
for shot_idx = 1:1
    no_discharge{shot_idx} = read_bgdata_shot(file_no_discharge, ...
        shot_idx, rotation_deg, gaussian_sigma);
end
no_discharge{2} = read_bgdata_shot(file_no_discharge_10psi, ...
    1, rotation_deg, gaussian_sigma);
no_discharge_with_filter{1} = read_bgdata_shot(...
    file_no_discharge_with_filter_7j, 1, rotation_deg, gaussian_sigma);
no_discharge_with_filter{2} = read_bgdata_shot(...
    file_no_discharge_with_filter, 1, rotation_deg, gaussian_sigma);

% First measure every centroid from the lower-left corner of the rotated
% image, then redefine the origin as the mean no-discharge centroid.
centroid_no_filter_abs = zeros(numel(no_filter), 2);
for shot_idx = 1:numel(no_filter)
    centroid_no_filter_abs(shot_idx, :) = centroid_from_lower_left(...
        no_filter{shot_idx}, centroid_cut.no_filter, centroid_min_area, ...
        centroid_locator_sigma, centroid_roi_half_width);
end
centroid_with_filter_abs = zeros(numel(with_filter), 2);
for shot_idx = 1:numel(with_filter)
    centroid_with_filter_abs(shot_idx, :) = ...
        centroid_from_lower_left(with_filter{shot_idx}, ...
        centroid_cut.with_filter, centroid_min_area, ...
        centroid_locator_sigma, centroid_roi_half_width);
end
centroid_no_discharge_abs = zeros(numel(no_discharge), 2);
for shot_idx = 1:numel(no_discharge)
    centroid_no_discharge_abs(shot_idx, :) = ...
        centroid_from_lower_left(no_discharge{shot_idx}, ...
        centroid_cut.no_discharge, centroid_min_area, ...
        centroid_locator_sigma, centroid_roi_half_width);
end
centroid_no_discharge_with_filter_abs = zeros(...
    numel(no_discharge_with_filter), 2);
for shot_idx = 1:numel(no_discharge_with_filter)
    centroid_no_discharge_with_filter_abs(shot_idx, :) = ...
        centroid_from_lower_left(no_discharge_with_filter{shot_idx}, ...
        centroid_cut.no_discharge_with_filter, centroid_min_area_weak, ...
        centroid_locator_sigma, centroid_roi_half_width);
end
% Define zero from the mean of the independently identified no-discharge
% spots. Averaging centroids avoids stationary detector texture becoming
% artificially prominent when the raw frames themselves are averaged.
reference_origin_um = mean(centroid_no_discharge_abs, 1);
[~, reference_shot_idx] = min(vecnorm(centroid_no_discharge_abs - ...
    reference_origin_um, 2, 2));
reference_image = no_discharge{reference_shot_idx};
reference_image_centroid = centroid_no_discharge_abs(reference_shot_idx, :) ...
    - reference_origin_um;

% All image panels share the maximum intensity of panel (a).
reference_intensity = max(no_filter{1}(:));

%% Figure
fg = figure(3);
set(fg, 'Color', 'w', 'Units', 'pixels', ...
    'Position', [200, 80, 1000, 820]);
layout = tiledlayout(fg, 2, 2, 'TileSpacing', 'compact', ...
    'Padding', 'compact');

base_map = jet(1000);
n_fade = 60;
fade_to_white = [linspace(base_map(n_fade, 1), 1, n_fade)', ...
    linspace(base_map(n_fade, 2), 1, n_fade)', ...
    linspace(base_map(n_fade, 3), 1, n_fade)'];
image_map = [flipud(fade_to_white); base_map];

%% (a) Full output without spectral filter
ax_a = nexttile(layout, 1);
plot_spatial_image(ax_a, no_filter{1}, reference_origin_um, ...
    reference_intensity, 1, image_map, fs);
plot_centroid_marker(ax_a, centroid_no_filter_abs(1, :) - reference_origin_um);
title(ax_a, 'Without filter');
add_panel_label(ax_a, '(a)', fs);

%% (b) MIR-dominated output with silicon filter
ax_b = nexttile(layout, 2);
plot_spatial_image(ax_b, with_filter{1}, reference_origin_um, ...
    reference_intensity, 10, image_map, fs);
plot_centroid_marker(ax_b, centroid_with_filter_abs(1, :) - reference_origin_um);
title(ax_b, 'With filter ($\times 10$)', 'Interpreter', 'latex');
add_panel_label(ax_b, '(b)', fs);
cb = colorbar(ax_b);
cb.Label.String = 'Normalized intensity';
cb.FontName = 'Times New Roman';
cb.FontSize = fs - 2;

%% (c) Reference without discharge
ax_c = nexttile(layout, 3);
plot_spatial_image(ax_c, reference_image, reference_origin_um, ...
    reference_intensity, 1, image_map, fs);
plot_centroid_marker(ax_c, reference_image_centroid);
title(ax_c, 'Without discharge');
add_panel_label(ax_c, '(c)', fs);

%% (d) Centroids calculated directly from the raw frames
ax_d = nexttile(layout, 4);
hold(ax_d, 'on');

centroid_no_filter = centroid_no_filter_abs - reference_origin_um;
centroid_with_filter = centroid_with_filter_abs - reference_origin_um;
centroid_no_discharge = centroid_no_discharge_abs - reference_origin_um;
centroid_no_discharge_with_filter = ...
    centroid_no_discharge_with_filter_abs - reference_origin_um;
fprintf('Reference origin [x,y] (um): %.3f, %.3f\n', reference_origin_um);
fprintf('Representative no-discharge shot in panel (c): %d\n', ...
    reference_shot_idx);
fprintf('Mean centered no-discharge centroid [x,y] (um): %.3f, %.3f\n', ...
    mean(centroid_no_discharge, 1));
disp('Centered no-discharge centroids [x,y] (um):');
disp(centroid_no_discharge);

p_reference = plot(ax_d, centroid_no_discharge(:, 1), ...
    centroid_no_discharge(:, 2), 'ko', 'MarkerFaceColor', 'k', ...
    'MarkerSize', 6, 'LineStyle', 'none');
p_no_filter = plot(ax_d, centroid_no_filter(:, 1), centroid_no_filter(:, 2), ...
    'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 7, 'LineStyle', 'none');
p_filter = plot(ax_d, centroid_with_filter(:, 1), ...
    centroid_with_filter(:, 2), 'ro', 'MarkerFaceColor', 'r', ...
    'MarkerSize', 7, 'LineStyle', 'none');
p_reference_filter = plot(ax_d, centroid_no_discharge_with_filter(:, 1), ...
    centroid_no_discharge_with_filter(:, 2), 'go', ...
    'MarkerFaceColor', 'g', 'MarkerSize', 7, 'LineStyle', 'none');

xlabel(ax_d, '$x~(\mu\mathrm{m})$', 'Interpreter', 'latex');
ylabel(ax_d, '$y~(\mu\mathrm{m})$', 'Interpreter', 'latex');
title(ax_d, 'Centroids from raw frames');
set(ax_d, 'FontName', 'Times New Roman', 'FontSize', fs, ...
    'Box', 'on', 'TickDir', 'out', 'TickLength', [0.02, 0.025]);
axis(ax_d, 'equal');
xlim(ax_d, [-300, 300]);
ylim(ax_d, [-300, 300]);
legend(ax_d, [p_reference, p_no_filter, p_filter, p_reference_filter], ...
    {'Without discharge', 'Without filter', 'With filter', ...
    'Without discharge + filter'}, ...
    'Location', 'southwest', 'Box', 'on', ...
    'FontName', 'Times New Roman', 'FontSize', fs - 3);
add_panel_label(ax_d, '(d)', fs);

%% Export
output_png = fullfile(output_dir, 'Fig3_4figs.png');
output_eps = fullfile(output_dir, 'Fig3_4figs.eps');
output_pdf = fullfile(output_dir, 'Fig3_4figs.pdf');
exportgraphics(fg, output_png, 'Resolution', 300);
print(fg, output_eps, '-depsc', '-painters');
exportgraphics(fg, output_pdf, 'ContentType', 'vector');
fprintf('Saved Fig. 3 to:\n%s\n%s\n%s\n', ...
    output_png, output_eps, output_pdf);

function image = read_bgdata_shot(file_path, shot_idx, rotation_deg, sigma)
    raw = h5read(file_path, sprintf('/BG_DATA/%d/DATA', shot_idx));
    height = double(h5read(file_path, ...
        sprintf('/BG_DATA/%d/RAWFRAME/HEIGHT', shot_idx)));
    width = double(h5read(file_path, ...
        sprintf('/BG_DATA/%d/RAWFRAME/WIDTH', shot_idx)));
    image = reshape(double(raw), width, height)' ./ 2^19;

    % Remove only the scalar detector offset; retain all spatial structure.
    image = max(image - median(image(:)), 0);
    image = interp2(image, 1, 'linear');
    image = imrotate(image, rotation_deg, 'bilinear', 'loose');
    image = imgaussfilt(image, sigma);
end

function plot_spatial_image(ax, image, reference_origin_um, ...
        reference_intensity, gain, map, fs)
    [x_axis, y_axis] = image_axes(image, reference_origin_um);
    display_image = gain * image / reference_intensity;
    imagesc(ax, x_axis, y_axis, display_image);
    set(ax, 'YDir', 'normal', 'FontName', 'Times New Roman', ...
        'FontSize', fs, 'Box', 'on', 'TickDir', 'out', ...
        'TickLength', [0.02, 0.025]);
    colormap(ax, map);
    caxis(ax, [0, 1]);
    axis(ax, 'equal');
    xlim(ax, [-300, 300]);
    ylim(ax, [-300, 300]);
    xlabel(ax, '$x~(\mu\mathrm{m})$', 'Interpreter', 'latex');
    ylabel(ax, '$y~(\mu\mathrm{m})$', 'Interpreter', 'latex');
end

function centroid = intensity_centroid(image, threshold_fraction, min_area, ...
        locator_sigma, roi_half_width)
    % Use a signal mask only for centroid extraction. The displayed image
    % remains unmasked so weak spatial components are retained.
    % Locate the broad spot automatically, independent of rotation angle.
    locator = imgaussfilt(image, locator_sigma);
    [~, locator_idx] = max(locator(:));
    [center_row, center_column] = ind2sub(size(image), locator_idx);
    half_rows = round(roi_half_width * size(image, 1));
    half_columns = round(roi_half_width * size(image, 2));
    roi = false(size(image));
    column_range = max(1, center_column - half_columns): ...
        min(size(image, 2), center_column + half_columns);
    row_range = max(1, center_row - half_rows): ...
        min(size(image, 1), center_row + half_rows);
    roi(row_range, column_range) = true;
    mask_source = image;
    mask_source(~roi) = 0;
    mask = mask_source >= threshold_fraction * max(mask_source(:));
    mask = bwareaopen(mask, min_area);
    components = bwconncomp(mask);
    if components.NumObjects > 1
        [~, peak_idx] = max(mask_source(:));
        strongest_idx = find(cellfun(@(idx) any(idx == peak_idx), ...
            components.PixelIdxList), 1);
        if isempty(strongest_idx)
            component_peak = cellfun(@(idx) max(mask_source(idx)), ...
                components.PixelIdxList);
            [~, strongest_idx] = max(component_peak);
        end
        mask(:) = false;
        mask(components.PixelIdxList{strongest_idx}) = true;
    end
    signal = image .* mask;
    [column_grid, row_grid] = meshgrid(1:size(image, 2), 1:size(image, 1));
    weight = sum(signal(:));
    if weight == 0
        error('No signal region survived the centroid threshold.');
    end
    centroid = [sum(row_grid(:) .* signal(:)) / weight, ...
        sum(column_grid(:) .* signal(:)) / weight];
end

function centroid_um = centroid_from_lower_left(image, ...
        threshold_fraction, min_area, locator_sigma, roi_half_width)
    centroid = intensity_centroid(image, threshold_fraction, min_area, ...
        locator_sigma, roi_half_width);
    scale_um = spatial_scale_um(image);
    centroid_um = [(centroid(2) - 1) * scale_um, ...
        (centroid(1) - 1) * scale_um];
end

function [x_axis, y_axis] = image_axes(image, reference_origin_um)
    scale_um = spatial_scale_um(image);
    x_axis = [0, size(image, 2) - 1] * scale_um - reference_origin_um(1);
    y_axis = [0, size(image, 1) - 1] * scale_um - reference_origin_um(2);
end

function scale_um = spatial_scale_um(image)
    % interp2(image,1) changes the 320-pixel detector axis to 639 samples.
    % Rotation with the 'loose' option adds padding but must not change the
    % physical pixel pitch.
    detector_samples_after_interp = 2 * 320 - 1;
    scale_um = 25.6 / detector_samples_after_interp / 4.5 * 250;
end

function add_panel_label(ax, label, fs)
    text(ax, -0.11, 1.015, label, 'Units', 'normalized', ...
        'FontName', 'Times New Roman', 'FontSize', fs + 4, ...
        'FontWeight', 'bold', 'VerticalAlignment', 'bottom', ...
        'HorizontalAlignment', 'left', 'Clipping', 'off');
end

function plot_centroid_marker(ax, centroid_um)
    hold(ax, 'on');
    plot(ax, centroid_um(1), centroid_um(2), 'w+', ...
        'MarkerSize', 12, 'LineWidth', 1.8);
    plot(ax, centroid_um(1), centroid_um(2), 'k+', ...
        'MarkerSize', 8, 'LineWidth', 1.0);
end

