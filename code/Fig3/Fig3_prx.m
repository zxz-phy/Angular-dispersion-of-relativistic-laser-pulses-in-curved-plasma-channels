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
file_with_filter = fullfile(data_dir, ...
    '20250718_OAP+透镜pyro成像_外部加硅片_34500fs18kV20psi_真红外信号.bgData');

%% Processing settings: no thresholding or connected-component masking
rotation_deg = -90;
gaussian_sigma = 1.0;
% Hard cuts are used only for centroid extraction. They do not alter the
% lightly processed images displayed in panels (a)-(c).
centroid_cut.no_filter = 0.80;
centroid_cut.with_filter = 0.80;
display_cut.with_filter = 0.30;
centroid_min_area = 50;
centroid_min_area_with_filter = 1;
centroid_locator_sigma = 12;
centroid_roi_half_width = 0.18;
fs = 20;
display_radius_um = 150;

% Lateral displacement through the tilted silicon plate. The plate normal
% forms angle_silicon with the incident beam; the correction is applied to
% the filtered image along +x before the image rotation.
d_silicon_um = 500;
n_silicon = 3.42;
angle_silicon_deg = 30;
angle_refraction_deg = asind(sind(angle_silicon_deg) / n_silicon);
silicon_lateral_shift_um = d_silicon_um * ...
    sind(angle_silicon_deg - angle_refraction_deg) / cosd(angle_refraction_deg);

no_filter = read_bgdata_shot(file_no_filter, 1, rotation_deg, gaussian_sigma, 0);
with_filter = read_bgdata_shot(file_with_filter, 1, rotation_deg, gaussian_sigma, ...
    silicon_lateral_shift_um);

% First measure every centroid from the lower-left corner of the rotated
% image, then redefine the origin as the mean no-discharge centroid.
centroid_no_filter_abs = centroid_from_lower_left(no_filter, ...
    centroid_cut.no_filter, centroid_min_area, ...
    centroid_locator_sigma, centroid_roi_half_width, true);
centroid_with_filter_abs = centroid_from_lower_left(with_filter, ...
    centroid_cut.with_filter, centroid_min_area_with_filter, ...
    centroid_locator_sigma, centroid_roi_half_width, true);
fit_no_filter_abs = fit_gaussian_spot(no_filter, centroid_locator_sigma, ...
    centroid_roi_half_width);
fit_with_filter_abs = fit_gaussian_spot(with_filter, centroid_locator_sigma, ...
    centroid_roi_half_width);
% Set the displayed unfiltered shot exactly at the spatial origin, then
% apply the same translation to every filtered and unfiltered shot.
reference_origin_um = centroid_no_filter_abs(1, :);

% Correct the silicon-filtered signal, then normalize both panels to the
% largest intensity across the unfiltered and corrected filtered images.
filter_attenuation_correction = 10;
reference_intensity = max([max(no_filter(:)), ...
    filter_attenuation_correction * max(with_filter(:))]);

%% Figure
fg = figure(3);
set(fg, 'Color', 'w', 'Units', 'pixels', ...
    'Position', [140, 140, 1200, 520]);
layout = tiledlayout(fg, 1, 2, 'TileSpacing', 'compact', ...
    'Padding', 'loose');

base_map = jet(1000);
n_fade = 60;
fade_to_white = [linspace(base_map(n_fade, 1), 1, n_fade)', ...
    linspace(base_map(n_fade, 2), 1, n_fade)', ...
    linspace(base_map(n_fade, 3), 1, n_fade)'];
image_map = [flipud(fade_to_white); base_map];

%% (a) Full output without spectral filter
ax_a = nexttile(layout, 1);
plot_spatial_image(ax_a, no_filter, reference_origin_um, ...
    reference_intensity, 1, image_map, fs, centroid_cut.no_filter, ...
    centroid_min_area, centroid_locator_sigma, centroid_roi_half_width, 0, ...
    display_radius_um);
plot_centroid_marker(ax_a, centroid_no_filter_abs - reference_origin_um);
plot_aperture_circle(ax_a, display_radius_um);
title(ax_a, 'Without filter');
add_panel_label(ax_a, '(a)', fs);

%% (b) MIR-dominated output with silicon filter
ax_b = nexttile(layout, 2);
plot_spatial_image(ax_b, with_filter, reference_origin_um, ...
    reference_intensity, filter_attenuation_correction, image_map, fs, centroid_cut.with_filter, ...
    centroid_min_area_with_filter, centroid_locator_sigma, ...
    centroid_roi_half_width, display_cut.with_filter, display_radius_um);
plot_centroid_marker(ax_b, centroid_with_filter_abs - reference_origin_um);
plot_aperture_circle(ax_b, display_radius_um);
title(ax_b, 'With filter ($\times 5$)', 'Interpreter', 'latex');
add_panel_label(ax_b, '(b)', fs);
cb = colorbar(ax_b);
cb.Label.String = 'Intensity (a. u.)';
cb.FontName = 'Times New Roman';
cb.FontSize = fs ;
cb.Ticks = [0, 0.25, 0.5, 0.75, 1];

fprintf('Reference origin [x,y] (um): %.3f, %.3f\n', reference_origin_um);
fprintf('Silicon plate: incidence = %.3f deg, refraction = %.3f deg, x shift = +%.3f um\n', ...
    angle_silicon_deg, angle_refraction_deg, silicon_lateral_shift_um);
fprintf('Panel (a) intensity-weighted centroid, absolute [x,y] (um): %.3f, %.3f\n', ...
    centroid_no_filter_abs);
fprintf('Filtered centroid relative to panel (a) [x,y] (um): %.3f, %.3f\n', ...
    centroid_with_filter_abs - reference_origin_um);
fprintf('Panel (b) intensity-weighted centroid, absolute [x,y] (um): %.3f, %.3f\n', ...
    centroid_with_filter_abs);
print_gaussian_fit('Panel (a)', fit_no_filter_abs, reference_origin_um);
print_gaussian_fit('Panel (b)', fit_with_filter_abs, reference_origin_um);

% Prevent MATLAB axes toolbars from appearing in exported publication figures.
ax_a.Toolbar.Visible = 'off';
ax_b.Toolbar.Visible = 'off';
drawnow;

%% Export
output_png = fullfile(output_dir, 'Fig3.png');
output_eps = fullfile(output_dir, 'Fig3.eps');
output_pdf = fullfile(output_dir, 'Fig3.pdf');
exportgraphics(fg, output_png, 'Resolution', 300);
print(fg, output_eps, '-depsc', '-painters');
exportgraphics(fg, output_pdf, 'ContentType', 'vector');
fprintf('Saved Fig. 3 to:\n%s\n%s\n%s\n', ...
    output_png, output_eps, output_pdf);

function image = read_bgdata_shot(file_path, shot_idx, rotation_deg, sigma, x_shift_um)
    raw = h5read(file_path, sprintf('/BG_DATA/%d/DATA', shot_idx));
    height = double(h5read(file_path, ...
        sprintf('/BG_DATA/%d/RAWFRAME/HEIGHT', shot_idx)));
    width = double(h5read(file_path, ...
        sprintf('/BG_DATA/%d/RAWFRAME/WIDTH', shot_idx)));
    image = reshape(double(raw), width, height)' ./ 2^19;

    % Remove only the scalar detector offset; retain all spatial structure.
    image = max(image - median(image(:)), 0);
    image = interp2(image, 1, 'linear');
    if x_shift_um ~= 0
        image = imtranslate(image, [x_shift_um / spatial_scale_um(image), 0], ...
            'OutputView', 'same', 'FillValues', 0);
    end
    % Restore the vertical detector orientation before rotating into the
    % climbing-mirror diagnostic coordinates.
    image = flipud(image);
    image = imrotate(image, rotation_deg, 'bilinear', 'loose');
    image = imgaussfilt(image, sigma);
end

function plot_spatial_image(ax, image, reference_origin_um, ...
        reference_intensity, signal_gain, map, fs, threshold_fraction, ...
        min_area, locator_sigma, roi_half_width, display_threshold, ...
        display_radius_um)
    [x_axis, y_axis] = image_axes(image, reference_origin_um);
    display_image = image / reference_intensity;
    if signal_gain > 1
        signal_mask = build_signal_mask(image, display_threshold, ...
            min_area, locator_sigma, roi_half_width);
        display_image(signal_mask) = signal_gain * display_image(signal_mask);
    end
    scale_um = spatial_scale_um(image);
    x_pixels = (0:size(image, 2) - 1) * scale_um - reference_origin_um(1);
    y_pixels = (0:size(image, 1) - 1) * scale_um - reference_origin_um(2);
    [x_grid, y_grid] = meshgrid(x_pixels, y_pixels);
    display_image(x_grid.^2 + y_grid.^2 > display_radius_um^2) = 0;
    imagesc(ax, x_axis, y_axis, display_image);
    set(ax, 'YDir', 'normal', 'FontName', 'Times New Roman', ...
        'FontSize', fs, 'Box', 'on', 'TickDir', 'out', ...
        'TickLength', [0.02, 0.025]);
    colormap(ax, map);
    caxis(ax, [0, 1]);
    axis(ax, 'equal');
    xlim(ax, [-150, 150]);
    ylim(ax, [-150, 150]);
    ax.XTick=[-150:50:150];
    ax.YTick=[-150:50:150];
    xlabel(ax, '$x~(\mu\mathrm{m})$', 'Interpreter', 'latex');
    ylabel(ax, '$y~(\mu\mathrm{m})$', 'Interpreter', 'latex');
end

function centroid = intensity_centroid(image, threshold_fraction, min_area, ...
        locator_sigma, roi_half_width, intensity_weighted)
    % Use a signal mask only for centroid extraction. The displayed image
    % remains unmasked so weak spatial components are retained.
    mask = build_signal_mask(image, threshold_fraction, min_area, ...
        locator_sigma, roi_half_width);
    if intensity_weighted
        signal = image .* mask;
    else
        signal = double(mask);
    end
    [column_grid, row_grid] = meshgrid(1:size(image, 2), 1:size(image, 1));
    weight = sum(signal(:));
    if weight == 0
        error('No signal region survived the centroid threshold.');
    end
    centroid = [sum(row_grid(:) .* signal(:)) / weight, ...
        sum(column_grid(:) .* signal(:)) / weight];
end

function mask = build_signal_mask(image, threshold_fraction, min_area, ...
        locator_sigma, roi_half_width)
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
    if ~any(mask(:))
        error('No signal region survived the centroid threshold.');
    end
end

function centroid_um = centroid_from_lower_left(image, ...
        threshold_fraction, min_area, locator_sigma, roi_half_width, intensity_weighted)
    centroid = intensity_centroid(image, threshold_fraction, min_area, ...
        locator_sigma, roi_half_width, intensity_weighted);
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
    % The two diameter measurements differ because the printed exit face is
    % not perfectly circular. Use their arithmetic mean as an isotropic scale.
    raw_pixel_scale_um = mean([300 / (228 - 159), 300 / (146 - 57)]);
    scale_um = raw_pixel_scale_um / 2;
end

function fit_result = fit_gaussian_spot(image, locator_sigma, roi_half_width)
    % Fit an elliptical 2D Gaussian within the spot-centered ROI.
    locator = imgaussfilt(image, locator_sigma);
    [~, locator_idx] = max(locator(:));
    [center_row, center_column] = ind2sub(size(image), locator_idx);
    half_rows = round(roi_half_width * size(image, 1));
    half_columns = round(roi_half_width * size(image, 2));
    rows = max(1, center_row-half_rows):min(size(image, 1), center_row+half_rows);
    columns = max(1, center_column-half_columns):min(size(image, 2), center_column+half_columns);
    fit_image = image(rows, columns);
    [column_grid, row_grid] = meshgrid(columns, rows);
    use = fit_image >= 0.10 * max(fit_image(:));
    x = column_grid(use);
    y = row_grid(use);
    z = fit_image(use);
    scale_um = spatial_scale_um(image);
    p0 = [max(z), center_column, center_row, log(12), log(12), 0];
    model = @(p) p(1) * exp(-0.5 * (((x-p(2))./exp(p(4))).^2 + ...
        ((y-p(3))./exp(p(5))).^2)) + p(6);
    options = optimset('Display', 'off', 'MaxFunEvals', 5000, 'MaxIter', 5000);
    p = fminsearch(@(p) sum((model(p)-z).^2), p0, options);
    sigma_x_um = exp(p(4)) * scale_um;
    sigma_y_um = exp(p(5)) * scale_um;
    fit_result.center_abs_um = [(p(2)-1)*scale_um, (p(3)-1)*scale_um];
    fit_result.w0_x_um = 2 * sigma_x_um;
    fit_result.w0_y_um = 2 * sigma_y_um;
    fit_result.fwhm_x_um = 2 * sqrt(2*log(2)) * sigma_x_um;
    fit_result.fwhm_y_um = 2 * sqrt(2*log(2)) * sigma_y_um;
end

function print_gaussian_fit(label, fit_result, reference_origin_um)
    fprintf('%s Gaussian-fit center, absolute [x,y] (um): %.3f, %.3f\n', ...
        label, fit_result.center_abs_um);
    fprintf('%s Gaussian-fit center, relative [x,y] (um): %.3f, %.3f\n', ...
        label, fit_result.center_abs_um-reference_origin_um);
    fprintf('%s Gaussian w0 [x,y] (um): %.3f, %.3f; FWHM [x,y] (um): %.3f, %.3f\n', ...
        label, fit_result.w0_x_um, fit_result.w0_y_um, ...
        fit_result.fwhm_x_um, fit_result.fwhm_y_um);
end

function add_panel_label(ax, label, fs)
    text(ax, -0.11, 1.045, label, 'Units', 'normalized', ...
        'FontName', 'Times New Roman', 'FontSize', fs + 6, ...
        'FontWeight', 'bold', 'VerticalAlignment', 'bottom', ...
        'HorizontalAlignment', 'left', 'Clipping', 'off');
end

function plot_aperture_circle(ax, radius_um)
    angle = linspace(0, 2 * pi, 500);
    hold(ax, 'on');
    plot(ax, radius_um * cos(angle), radius_um * sin(angle), ...
        'k-', 'LineWidth', 0.8);
end

function plot_centroid_marker(ax, centroid_um)
    hold(ax, 'on');
    plot(ax, centroid_um(1), centroid_um(2), 'w+', ...
        'MarkerSize', 12, 'LineWidth', 1.8);
    plot(ax, centroid_um(1), centroid_um(2), 'k+', ...
        'MarkerSize', 8, 'LineWidth', 1.0);
end
