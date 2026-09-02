clc;
clear;
close all;

%% User controls
% Vertex order: [top-left, top-right, bottom-left, bottom-right].
rainbow_x = [0, 0.01, 0.1, 0.4];
rainbow_y = [0, 0, -1, -1];
rainbow_common_center = [0, 1];
len = 1e2;
white_mode = 'on'; % 'on' hides the axes box, ticks, and labels.

%% Corresponding points on the top and bottom edges
top_x = linspace(rainbow_x(1), rainbow_x(2), len);
top_y = linspace(rainbow_y(1), rainbow_y(2), len);
bottom_x = linspace(rainbow_x(3), rainbow_x(4), len);
bottom_y = linspace(rainbow_y(3), rainbow_y(4), len);

%% Figure and rainbow colors
fg = figure(1);
set(fg, 'Color', 'w', 'Units', 'pixels', 'Position', [300, 160, 800, 600]);
ax = axes(fg);
hold(ax, 'on');

% Blue -> cyan -> yellow -> red, matching the wavelength colors in Fig. 2(a).
rainbow_map = wavelength_colormap(len);
n_arc_points = 400;

for line_idx = 1:len
    start_point = [top_x(line_idx), top_y(line_idx)];
    end_point = [bottom_x(line_idx), bottom_y(line_idx)];

    start_vector = start_point - rainbow_common_center;
    end_vector = end_point - rainbow_common_center;
    start_radius = hypot(start_vector(1), start_vector(2));
    end_radius = hypot(end_vector(1), end_vector(2));
    start_angle = atan2(start_vector(2), start_vector(1));
    end_angle = atan2(end_vector(2), end_vector(1));

    % Select the shorter continuous angular path. The radius interpolation
    % makes the curve pass through both endpoints. If the endpoint radii are
    % equal, this is an exact circular arc about rainbow_common_center.
    angle_pair = unwrap([start_angle, end_angle]);
    curve_parameter = linspace(0, 1, n_arc_points);
    curve_angle = angle_pair(1) + ...
        curve_parameter * (angle_pair(2) - angle_pair(1));
    curve_radius = start_radius + ...
        curve_parameter * (end_radius - start_radius);

    curve_x = rainbow_common_center(1) + curve_radius .* cos(curve_angle);
    curve_y = rainbow_common_center(2) + curve_radius .* sin(curve_angle);
    plot(ax, curve_x, curve_y, 'Color', rainbow_map(line_idx, :), ...
        'LineWidth', 2.5);
end

%% Axes appearance
all_x = [rainbow_x, rainbow_common_center(1)];
all_y = [rainbow_y, rainbow_common_center(2)];
x_padding = 0.08 * max(range(all_x), eps);
y_padding = 0.08 * max(range(all_y), eps);
xlim(ax, [min(rainbow_x)-x_padding, max(rainbow_x)+x_padding]);
ylim(ax, [min(rainbow_y)-y_padding, max(rainbow_y)+y_padding]);
axis(ax, 'equal');
set(ax, 'FontName', 'Times New Roman', 'FontSize', 16, ...
    'Box', 'on', 'TickDir', 'out');
xlabel(ax, 'x');
ylabel(ax, 'y');

if strcmpi(white_mode, 'on')
    axis(ax, 'off');
end

function map = wavelength_colormap(n)
anchors = [0.00, 0.00, 0.75; ...
           0.00, 0.65, 1.00; ...
           0.95, 0.90, 0.15; ...
           1.00, 0.30, 0.00; ...
           0.70, 0.00, 0.00];
anchor_x = linspace(0, 1, size(anchors, 1));
map = interp1(anchor_x, anchors, linspace(0, 1, n), 'pchip');
map = min(max(map, 0), 1);
end
