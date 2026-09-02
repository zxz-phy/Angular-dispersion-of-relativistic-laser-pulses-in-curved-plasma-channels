clc; clear; close all;

script_dir = fileparts(mfilename('fullpath'));
figs_root = fileparts(fileparts(script_dir));
data_dir = fullfile(figs_root, 'data', 'Fig4');
output_dir = fullfile(figs_root, 'Figs', 'Fig4');
if ~isfolder(output_dir), mkdir(output_dir); end

red = load(fullfile(data_dir, ...
    'Prism_Thick_dispersion_relation_L3mm_S2000R20000_1e3lines.mat'));
blue = load(fullfile(data_dir, ...
    'Prism_Thick_dispersion_relation_L3mm_S6000R20000_1e3lines.mat'));
trajectory = load(fullfile(data_dir, ...
    'Prism_Thick_dispersion_relation_L3mm_S10000R20000_all_data.mat'), ...
    'x1new', 'y1new', 'cneed', 'lam_prism', 'L_straight', 'S', 'R');
lam_red = red.lam_prism;
lam_blue = blue.lam_prism;
% Convert laboratory-frame angles to deviations relative to the channel axis:
% theta = theta_lab - theta_0, with theta_0 = S/R.
theta0_red = rad2deg(2 / 20);   % S = 2 mm, R = 20 mm
theta0_blue = rad2deg(6 / 20);  % S = 6 mm, R = 20 mm
theta_red = -red.angle_prism - theta0_red;
theta_blue = -blue.angle_prism - theta0_blue;

% Apply the same physical-width Savitzky-Golay smoothing to both cases,
% then use a centered numerical derivative. This avoids amplifying the
% small point-to-point structure differently in panels (b) and (c).
derivative_smoothing_window_um = 1.0;
derivative_post_smoothing_window_um = 2.5;
[theta_red_s, dtheta_red] = smooth_and_differentiate( ...
    lam_red, theta_red, derivative_smoothing_window_um, ...
    derivative_post_smoothing_window_um);
[theta_blue_s, dtheta_blue] = smooth_and_differentiate( ...
    lam_blue, theta_blue, derivative_smoothing_window_um, ...
    derivative_post_smoothing_window_um);
dlam_red = lam_red;
dlam_blue = lam_blue;

fs = 18;
lambda_cut_a = [0.8, 10]; % um; wavelength-colorbar range in panel (a)
rotation_deg0_position = 8; % mm; channel tangent here is horizontal after rotation
ax_acut_1 = [7, 8, -0.6, -0.4]; % mm: [z_min z_max y_min y_max]
ax_acut_1_position = [0.3, 0.62, 0.2, 0.15]; % normalized figure position
ax_acut_2 = [11.2, 12.2, -1.05, -0.85]; % mm: [z_min z_max y_min y_max]
ax_acut_2_position = [0.55, 0.62, 0.2, 0.15]; % normalized figure position
fg = figure(4);
set(fg, 'Color', 'w', 'Units', 'pixels', 'Position', [260, 40, 760, 980]);

%% (a) Wavelength-dependent trajectories in the S = 10 mm channel
ax_a = axes(fg, 'Position', [0.14, 0.685, 0.624, 0.255]);

% At z0 in the curved section, the centerline tangent angle is
% -(z0-L)/R. Rotate the plotted content by the opposite angle about that
% centerline point; the axes box itself remains horizontal and unrotated.
L_mm = trajectory.L_straight / 1e3;
S_mm = trajectory.S / 1e3;
R_mm = trajectory.R / 1e3;
if rotation_deg0_position < L_mm || rotation_deg0_position > L_mm + S_mm
    error('rotation_deg0_position must lie within the channel.');
end
arc_at_rotation = max(rotation_deg0_position - L_mm, 0);
rotation_angle = arc_at_rotation / R_mm;
if rotation_deg0_position <= L_mm
    pivot_x = rotation_deg0_position;
    pivot_y = 0;
else
    pivot_x = L_mm + R_mm * sin(rotation_angle);
    pivot_y = -R_mm + R_mm * cos(rotation_angle);
end

lambda_map = wavelength_colormap(256);
lambda_min = lambda_cut_a(1);
lambda_max = lambda_cut_a(2);
% Apply lambda_cut_a consistently to both trajectories and colorbar.
selected_idx = find(trajectory.lam_prism >= lambda_min & ...
    trajectory.lam_prism <= lambda_max);
if isempty(selected_idx)
    error('No panel-(a) trajectories fall inside lambda_cut_a.');
end
[~, wavelength_order] = sort(trajectory.lam_prism(selected_idx), 'descend');
selected_idx = selected_idx(wavelength_order);

% Rotate the selected trajectories first, then derive the exact horizontal
% cut from their rotated start-to-end extent.
rotated_x = cell(size(selected_idx));
rotated_y = cell(size(selected_idx));
all_rotated_x = [];
all_rotated_y = [];
for selected_number = 1:numel(selected_idx)
    line_idx = selected_idx(selected_number);
    [rotated_x{selected_number}, rotated_y{selected_number}] = ...
        rotate_about_point(trajectory.x1new(line_idx,:)/1e3, ...
        trajectory.y1new(line_idx,:)/1e3, pivot_x, pivot_y, rotation_angle);
    all_rotated_x = [all_rotated_x, rotated_x{selected_number}]; %#ok<AGROW>
    all_rotated_y = [all_rotated_y, rotated_y{selected_number}]; %#ok<AGROW>
end
z_limits = [min(all_rotated_x), max(all_rotated_x)];
y_padding = 0.55;
y_limits_a = [min(all_rotated_y)-y_padding, max(all_rotated_y)+y_padding];

% Build the density directly on a regular grid in the final rotated frame.
% Inverse-map each display pixel to the original straight-plus-arc geometry;
% this follows figure(23)'s density expression without rotating a coarse grid.
z_display = linspace(z_limits(1), z_limits(2), 1000);
y_display = linspace(y_limits_a(1), y_limits_a(2), 360);
[z_grid, y_grid] = meshgrid(z_display, y_display);
[original_x, original_y] = rotate_about_point(z_grid, y_grid, ...
    pivot_x, pivot_y, -rotation_angle);
density_offset_um = zeros(size(original_x));
valid_density = false(size(original_x));
straight_mask = original_x >= 0 & original_x <= L_mm;
density_offset_um(straight_mask) = original_y(straight_mask) * 1e3;
valid_density(straight_mask) = true;
curve_angle = atan2(original_x-L_mm, original_y+R_mm);
curved_mask = curve_angle >= 0 & curve_angle <= S_mm/R_mm;
density_offset_um(curved_mask) = ...
    (sqrt((original_x(curved_mask)-L_mm).^2 + ...
    (original_y(curved_mask)+R_mm).^2)-R_mm) * 1e3;
valid_density(curved_mask) = true;
ne = 1 + 4*density_offset_um.^2/16^4;
ne(ne >= 15 | ~valid_density) = 0;
imagesc(ax_a, z_display, y_display, ne);
set(ax_a, 'YDir', 'normal'); hold(ax_a, 'on');

% Draw from long to short wavelength, so the blue short-wavelength
% trajectory remains visible on top where the paths overlap.
for selected_number = 1:numel(selected_idx)
    line_idx = selected_idx(selected_number);
    lambda_norm = (trajectory.lam_prism(line_idx)-lambda_min) / ...
        (lambda_max-lambda_min);
    lambda_norm = min(max(lambda_norm,0),1);
    color_idx = 1 + round(lambda_norm*(size(lambda_map,1)-1));
    plot(ax_a, rotated_x{selected_number}, rotated_y{selected_number}, ...
        'Color', lambda_map(color_idx,:), 'LineWidth', 1.15);
end
xlim(ax_a, z_limits);
ylim(ax_a, y_limits_a);
colormap(ax_a, flipud(gray(64))); caxis(ax_a, [0, max(ne(:))]);

% Draw the z scale directly in panel-(a) data coordinates. With unequal
% display scales sz and sy, its on-screen angle automatically becomes
% beta = atan((sy/sz)*tan(rotation_angle)), exactly like the straight section.
z_span = diff(z_limits);
y_span = diff(y_limits_a);
axis(ax_a, 'off');
scale_origin_x = z_limits(1) + 0.075*z_span;
scale_origin_y = y_limits_a(2) - 0.22*y_span;
scale_z = 1.0; % mm
z_vector = scale_z*[cos(rotation_angle), sin(rotation_angle)];
quiver(ax_a, scale_origin_x, scale_origin_y, z_vector(1), z_vector(2), 0, ...
    'Color', 'k', 'LineWidth', 2.0, 'MaxHeadSize', 0.35);
pixels_per_z = fg.Position(3)*ax_a.Position(3)/z_span;
pixels_per_y = fg.Position(4)*ax_a.Position(4)/y_span;
z_display_angle_deg = rad2deg(atan2( ...
    sin(rotation_angle)*pixels_per_y, ...
    cos(rotation_angle)*pixels_per_z));
text(ax_a, scale_origin_x+0.5*z_vector(1), ...
    scale_origin_y+0.5*z_vector(2)-0.10, ...
    '1 mm', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
    'Rotation', z_display_angle_deg, ...
    'FontName', 'Times New Roman', 'FontSize', fs-2);
text(ax_a, scale_origin_x+z_vector(1)+0.08, ...
    scale_origin_y+z_vector(2), '$z$', ...
    'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', fs, ...
    'Rotation', z_display_angle_deg);

% Use a single compact horizontal wavelength colorbar at the lower left.
ax_lambda_cb = axes(fg, 'Position', [0.14, 0.685, 0.624, 0.255], ...
    'Visible', 'off', 'HitTest', 'off');
colormap(ax_lambda_cb, lambda_map);
caxis(ax_lambda_cb, [lambda_min, lambda_max]);
cb_lambda = colorbar(ax_lambda_cb, 'southoutside');
cb_lambda.Position = [0.155, 0.670, 0.125, 0.014];
cb_lambda.Ticks = [0.8, 5, 10];
cb_lambda.Label.String = '\lambda (\mum)';
cb_lambda.Label.Interpreter = 'tex'; cb_lambda.Label.FontSize = fs;
add_panel_label_top(ax_a, '(a)', fs);

%% Two magnified trajectory windows linked to their panel-(a) regions
zoom_axes = gobjects(1,2);
zoom_boxes = {ax_acut_1, ax_acut_2};
zoom_positions = {ax_acut_1_position, ax_acut_2_position};
for zoom_idx = 1:2
    zoom_box = zoom_boxes{zoom_idx};
    zoom_position = zoom_positions{zoom_idx};

    rectangle(ax_a, 'Position', [zoom_box(1), zoom_box(3), ...
        zoom_box(2)-zoom_box(1), zoom_box(4)-zoom_box(3)], ...
        'EdgeColor', 'k', 'LineWidth', 1.4, 'LineStyle', '-');

    zoom_axes(zoom_idx) = axes(fg, 'Position', zoom_position);
    imagesc(zoom_axes(zoom_idx), z_display, y_display, ne);
    set(zoom_axes(zoom_idx), 'YDir', 'normal');
    hold(zoom_axes(zoom_idx), 'on');
    for selected_number = 1:numel(selected_idx)
        line_idx = selected_idx(selected_number);
        lambda_norm = (trajectory.lam_prism(line_idx)-lambda_min) / ...
            (lambda_max-lambda_min);
        lambda_norm = min(max(lambda_norm,0),1);
        color_idx = 1 + round(lambda_norm*(size(lambda_map,1)-1));
        plot(zoom_axes(zoom_idx), rotated_x{selected_number}, ...
            rotated_y{selected_number}, 'Color', lambda_map(color_idx,:), ...
            'LineWidth', 1.0);
    end
    xlim(zoom_axes(zoom_idx), zoom_box(1:2));
    ylim(zoom_axes(zoom_idx), zoom_box(3:4));
    colormap(zoom_axes(zoom_idx), flipud(gray(64)));
    caxis(zoom_axes(zoom_idx), [0, max(ne(:))]);
    set(zoom_axes(zoom_idx), 'FontName', 'Times New Roman', ...
        'FontSize', fs-7, 'Box', 'on', 'TickDir', 'out', ...
        'Layer', 'top', 'XTick', [], 'YTick', []);

    % Connect both lower corners of the ROI to the corresponding upper
    % corners of the magnified axes.
    [roi_left_x, roi_bottom_y] = data_to_figure(ax_a, ...
        zoom_box(1), zoom_box(3));
    [roi_right_x, ~] = data_to_figure(ax_a, ...
        zoom_box(2), zoom_box(3));
    annotation(fg, 'line', [roi_left_x, zoom_position(1)], ...
        [roi_bottom_y, zoom_position(2)+zoom_position(4)], ...
        'Color', 'k', 'LineWidth', 1.0);
    annotation(fg, 'line', [roi_right_x, zoom_position(1)+zoom_position(3)], ...
        [roi_bottom_y, zoom_position(2)+zoom_position(4)], ...
        'Color', 'k', 'LineWidth', 1.0);
end

%% Compact, vertically aligned line-plot group: panels (b)-(d)
ax_b = axes(fg, 'Position', [0.14, 0.455, 0.624, 0.135]);
plot(ax_b, lam_red, theta_red_s, 'r-', 'LineWidth', 3);
xlim(ax_b, [0.8, 10]); ylim(ax_b, [-3.9, -2.4]); yticks(ax_b, [-3.5, -3, -2.5]);
ylabel(ax_b, '$\theta~(^{\circ})$', 'Interpreter', 'latex');
legend(ax_b, '$S=2$ mm, $R=20$ mm', 'Interpreter', 'latex', ...
    'Location', 'southeast', 'Box', 'off');
add_panel_label(ax_b, '(b)', fs);

ax_c = axes(fg, 'Position', [0.14, 0.320, 0.624, 0.135]);
plot(ax_c, lam_blue, theta_blue_s, 'b-', 'LineWidth', 3);
xlim(ax_c, [0.8, 10]); ylim(ax_c, [1.9, 3.4]); yticks(ax_c, [1.5, 2, 2.5, 3]);
ylabel(ax_c, '$\theta~(^{\circ})$', 'Interpreter', 'latex');
legend(ax_c, '$S=6$ mm, $R=20$ mm', 'Interpreter', 'latex', ...
    'Location', 'northeast', 'Box', 'off');
add_panel_label(ax_c, '(c)', fs);

ax_d = axes(fg, 'Position', [0.14, 0.185, 0.624, 0.135]); hold(ax_d, 'on');
plot(ax_d, dlam_red, dtheta_red, 'r-', 'LineWidth', 3);
plot(ax_d, dlam_blue, dtheta_blue, 'b-', 'LineWidth', 3);
xlim(ax_d, [0.8, 10]); ylim(ax_d, [-0.6, 0.6]); yticks(ax_d, [-0.4, 0, 0.4]);
% xticks(ax_d, [0.8 2:2:10]);
xlabel(ax_d, '$\lambda~(\mu\mathrm{m})$', 'Interpreter', 'latex');
ylabel(ax_d, '$d\theta/d\lambda$', 'Interpreter', 'latex');
add_panel_label(ax_d, '(d)', fs);

set([ax_a, ax_b, ax_c, ax_d], 'FontName', 'Times New Roman', 'FontSize', fs, ...
    'Box', 'on', 'TickDir', 'out');
set([ax_b, ax_c, ax_d], ...
    'Box', 'on', 'TickDir', 'out', 'XTick', [0.8 2:2:10]);
set([ax_b, ax_c], 'XTickLabel', []);
set([ax_b.YLabel, ax_c.YLabel, ax_d.YLabel], 'Units', 'normalized');
ax_b.YLabel.Position = [-0.105, 0.5, 0];
ax_c.YLabel.Position = [-0.105, 0.5, 0];
ax_d.YLabel.Position = [-0.105, 0.5, 0];
uistack(zoom_axes, 'top');

% Report reproducible values used to describe panels (b)--(d).
lambda_report_um = [0.8, 2.0, 3.0, 10.0];
theta_red_report = interp1(lam_red, theta_red_s, lambda_report_um, 'pchip');
theta_blue_report = interp1(lam_blue, theta_blue_s, lambda_report_um, 'pchip');
dtheta_red_3um = interp1(dlam_red, dtheta_red, 3.0, 'pchip');
dtheta_blue_3um = interp1(dlam_blue, dtheta_blue, 3.0, 'pchip');
fprintf('Fig. 4 relative deviation angles theta (deg)\n');
fprintf('  lambda (um):  %7.3f %7.3f %7.3f %7.3f\n', lambda_report_um);
fprintf('  S = 2 mm:     %7.4f %7.4f %7.4f %7.4f\n', theta_red_report);
fprintf('  S = 6 mm:     %7.4f %7.4f %7.4f %7.4f\n', theta_blue_report);
fprintf('  S = 2 mm, theta(3.0)-theta(0.8) = %+0.4f deg\n', ...
    theta_red_report(3) - theta_red_report(1));
fprintf('  dtheta/dlambda at 3.0 um: S = 2 mm, %+0.4f deg/um; ', ...
    dtheta_red_3um);
fprintf('S = 6 mm, %+0.4f deg/um\n', dtheta_blue_3um);

exportgraphics(fg, fullfile(output_dir, 'Fig4.png'), 'Resolution', 600);
print(fg, fullfile(output_dir, 'Fig4.eps'), '-depsc', '-opengl', '-r600');
exportgraphics(fg, fullfile(output_dir, 'Fig4.pdf'), 'ContentType', 'vector');
fprintf('Saved Fig. 4 to %s\n', output_dir);

function add_panel_label(ax, label, fs)
text(ax, 0.015, 0.20, label, 'Units', 'normalized', ...
    'FontName', 'Times New Roman', 'FontSize', fs + 4, ...
    'FontWeight', 'bold', 'Clipping', 'off');
end

function add_panel_label_top(ax, label, fs)
text(ax, -0.035, 1.035, label, 'Units', 'normalized', ...
    'FontName', 'Times New Roman', 'FontSize', fs + 4, ...
    'FontWeight', 'bold', 'VerticalAlignment', 'bottom', 'Clipping', 'off');
end

function [x_rot, y_rot] = rotate_about_point(x, y, pivot_x, pivot_y, angle)
x_shift = x - pivot_x;
y_shift = y - pivot_y;
x_rot = pivot_x + x_shift * cos(angle) - y_shift * sin(angle);
y_rot = pivot_y + x_shift * sin(angle) + y_shift * cos(angle);
end

function map = wavelength_colormap(n)
% Balanced blue-cyan-yellow-red map with a substantial warm-color range.
anchors = [0.00, 0.00, 0.75; ...
           0.00, 0.65, 1.00; ...
           0.95, 0.90, 0.15; ...
           1.00, 0.30, 0.00; ...
           0.70, 0.00, 0.00];
anchor_x = linspace(0,1,size(anchors,1));
map = interp1(anchor_x, anchors, linspace(0,1,n), 'pchip');
map = min(max(map,0),1);
end

function [theta_smooth, derivative] = smooth_and_differentiate( ...
    lambda, theta, window_um, derivative_window_um)
lambda_column = lambda(:);
theta_column = theta(:);
delta_lambda = median(diff(lambda_column));
span = round(window_um/delta_lambda);
span = max(span, 5);
if mod(span,2) == 0
    span = span + 1;
end
span = min(span, numel(theta_column) - 1 + mod(numel(theta_column),2));
theta_smooth_column = smooth(theta_column, span, 'sgolay', 3);
derivative_column = gradient(theta_smooth_column, lambda_column);
derivative_span = round(derivative_window_um/delta_lambda);
derivative_span = max(derivative_span, 5);
if mod(derivative_span,2) == 0
    derivative_span = derivative_span + 1;
end
derivative_span = min(derivative_span, ...
    numel(theta_column) - 1 + mod(numel(theta_column),2));
derivative_column = smooth(derivative_column, derivative_span, 'sgolay', 2);
theta_smooth = reshape(theta_smooth_column, size(theta));
derivative = reshape(derivative_column, size(theta));
end

function [figure_x, figure_y] = data_to_figure(ax, data_x, data_y)
axes_position = ax.Position;
x_limits = ax.XLim;
y_limits = ax.YLim;
figure_x = axes_position(1) + ...
    (data_x-x_limits(1))/diff(x_limits)*axes_position(3);
figure_y = axes_position(2) + ...
    (data_y-y_limits(1))/diff(y_limits)*axes_position(4);
end
