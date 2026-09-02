clc; clear; close all;

if_export = false;

script_dir = fileparts(mfilename('fullpath'));
figs_root = fileparts(fileparts(script_dir));
data_dir = fullfile(figs_root, 'data', 'Fig5');
output_dir = fullfile(figs_root, 'Figs', 'Fig5');
ax_a_position=[0.250, 0.580,  0.50, 0.360];
ax_b_position=[0.250, 0.100,  0.50, 0.360];
cb_a_Position=[0.800, 0.580, 0.025, 0.360];
cb_b_Position=[0.800, 0.100, 0.025, 0.360];

if ~isfolder(output_dir), mkdir(output_dir); end
data = load(fullfile(data_dir, 'Thick_prism_heatmap_combined_20260821_0o8to0.mat'));

fg = figure(5);
set(fg, 'Color', 'w', 'Units', 'pixels', 'Position', [180, 30, 760, 1180]);
fs = 26; n = 100;
map = [linspace(0,1,n)', linspace(0,1,n)', ones(n,1); ...
       ones(n,1), linspace(1,0,n)', linspace(1,0,n)'];

ax_a = axes(fg, 'Position', ax_a_position);
imagesc(ax_a, data.S_query/1e3, data.R_query/1e3, data.theta_relative_need_heat_map_deg);
format_map_axes(ax_a, fs); colormap(ax_a, map); caxis(ax_a, [-10, 10]);
cb_a = colorbar(ax_a); cb_a.Position = cb_a_Position;
cb_a.Label.String = '\theta (°)';
cb_a.Label.Interpreter = 'tex'; cb_a.Label.FontName = 'Times New Roman';
cb_a.Label.FontSize = 20;
ax_a.Position = ax_a_position;
% pbaspect(ax_a, [1, 1, 1]);
hold(ax_a, 'on'); mask = data.S_query/1e3 >= 0.5;
contour(ax_a, data.S_query(mask)/1e3, data.R_query/1e3, ...
    data.theta_relative_need_heat_map_deg(:,mask), [0 0], 'k--', 'LineWidth', 2);
panel(ax_a, '(a)', fs);

ax_b = axes(fg, 'Position', ax_b_position);
dispersion_map = -data.dtdl_heat_map;
imagesc(ax_b, data.S_query/1e3, data.R_query/1e3, dispersion_map);
format_map_axes(ax_b, fs); colormap(ax_b, map); caxis(ax_b, [-1, 1]);
cb_b = colorbar(ax_b); cb_b.Position = cb_b_Position;
cb_b.Label.String = 'd\theta/d\lambda (°/\mum)';
cb_b.Label.Interpreter = 'tex'; cb_b.Label.FontName = 'Times New Roman';
cb_b.Label.FontSize = 20;
ax_b.Position = ax_b_position;
% pbaspect(ax_b, [1, 1, 1]);
hold(ax_b, 'on');
% Retain the original left transition, for which the weak-dispersion data
% do not define a visually stable zero contour.
x_left = linspace(1.32, 2.47, 1e4);
y_left = 40 * (x_left - 0.9).^6 + 9.5;
plot(ax_b, x_left, y_left, 'k--', 'LineWidth', 2);

% Use the zero contours calculated from the updated map for the two
% well-resolved transitions on the right.
plot_right_zero_contours(ax_b, data.S_query/1e3, data.R_query/1e3, ...
    dispersion_map, 3.0);

text(ax_b, 3.45, 44.5, {'Positive', 'dispersion'}, ...
    'HorizontalAlignment', 'center', 'FontName', 'Times New Roman', ...
    'FontSize', fs-8);
text(ax_b, 7.05, 44.5, {'Negative', 'dispersion'}, ...
    'HorizontalAlignment', 'center', 'FontName', 'Times New Roman', ...
    'FontSize', fs-8);
panel(ax_b, '(b)', fs);

if if_export
    exportgraphics(fg, fullfile(output_dir,'Fig5.png'), 'Resolution', 600);
    print(fg, fullfile(output_dir,'Fig5.eps'), '-depsc', '-painters');
    exportgraphics(fg, fullfile(output_dir,'Fig5.pdf'), 'ContentType', 'vector');
    fprintf('Saved Fig. 5 to %s\n', output_dir);
else
    fprintf('Fig. 5 export skipped (if_export = false).\n');
end

function format_map_axes(ax, fs)
set(ax,'YDir','normal','FontName','Times New Roman','FontSize',fs, ...
    'XLim',[0 10],'YLim',[10 50],'XTick',0:2:10,'YTick',10:10:50);
xlabel(ax,'S (mm)'); ylabel(ax,'R (mm)');
end
function panel(ax,label,fs)
text(ax,.025,.985,label,'Units','normalized','FontName','Times New Roman', ...
    'FontSize',fs,'FontWeight','bold','VerticalAlignment','top');
end

function plot_right_zero_contours(ax, x, y, map_data, minimum_s)
contour_matrix = contourc(x, y, map_data, [0 0]);
column = 1;
while column < size(contour_matrix, 2)
    point_count = contour_matrix(2, column);
    points = contour_matrix(:, column + (1:point_count));
    if median(points(1, :)) >= minimum_s
        plot(ax, points(1, :), points(2, :), 'k--', 'LineWidth', 2);
    end
    column = column + point_count + 1;
end
end
