clc;
clear;
close all;

%% User controls
t = 5200;                    % OSIRIS output time used in the main text (um)
if_2d_visual = true;         % Figure 2: source 2D field with disk indices
number_of_disks = 18;
disk_visualized = [7:2:14 15:18];
% Initial 5200-um inspection values, shifted with the moving simulation
% window relative to the 3100-um frame. These will be refined later.
disk_split_um = 5119.45;
leading_search_limit_um = 5129.78;
merged_lobe_center_um = 5117.55;
merged_lobe_range_um = [5115.78 5119.45];
trailing_search_limit_um = 5114.78;

% Display window around the laser centroid. Empty values trigger automatic
% limits around the intensity-weighted centroid.
x_window_um = 70;
y_window_um = 120;           % retain the complete deflected pulse, without blank domain

% Pseudo-3D display parameters. The third coordinate is a visualization
% coordinate and must not be interpreted as a simulated spatial dimension.
laser_ne_offset = 5;         % laser center above the ne=n0 plane
ne_amplitude = 10;           % vertical gain of (ne-n0) in Figure 1
laser_plane_z = laser_ne_offset;
laser_disk_z_radius = 3.4;
maximum_disk_half_width = 0.75;
disk_longitudinal_fatness = 2; % multiplier along propagation direction x
ne_alpha = 0.90;

% Figure-1 view angles [azimuth elevation] in degrees. Axis limits are
% automatic in this first 5200-um inspection.
fig1_view = [0 90];
surface_stride = 4;

%% Paths
code_dir = fileparts(mfilename('fullpath'));
figs_dir = fileparts(fileparts(code_dir));
data_dir = fullfile(figs_dir, 'data', 'Fig1');
e3_file = fullfile(data_dir, 'e3', ...
    sprintf('e3-%06d.h5', t));
charge_file = fullfile(data_dir, 'charge', ...
    sprintf('charge-Plasma_elec-%06d.h5', t));

assert(isfile(e3_file), 'Missing laser-field file: %s', e3_file);
assert(isfile(charge_file), 'Missing density file: %s', charge_file);

%% Read the 2D PIC data
% OSIRIS stores arrays as [x,y]. Transposition gives matrices indexed as
% [y,x], consistent with imagesc(x,y,data).
e3 = double(h5read(e3_file, '/e3'))' ./ 7.854;
ne = -double(h5read(charge_file, '/charge'))' .* 5.658;

xmin = double(h5readatt(e3_file, '/', 'XMIN'));
xmax = double(h5readatt(e3_file, '/', 'XMAX'));
x = linspace(xmin(1), xmax(1), size(e3, 2));
y = linspace(xmin(2), xmax(2), size(e3, 1));

% Use positive electron density for visualization and suppress isolated
% numerical undershoots.
ne = max(ne, 0);

%% Locate the pulse and crop the common field-density window
intensity = e3.^2;
longitudinal_energy = sum(intensity, 1);
[~, ix_centroid] = max(longitudinal_energy);
x_centroid = x(ix_centroid);

transverse_energy = sum(intensity, 2);
y_centroid = sum(y(:) .* transverse_energy) / max(sum(transverse_energy), eps);

x_limits = [x_centroid - 0.72*x_window_um, ...
            x_centroid + 0.28*x_window_um];
y_limits = [y_centroid - 0.5*y_window_um, ...
            y_centroid + 0.5*y_window_um];
x_limits(1) = max(x_limits(1), x(1));
x_limits(2) = min(x_limits(2), x(end));
y_limits(1) = max(y_limits(1), y(1));
y_limits(2) = min(y_limits(2), y(end));

ix = find(x >= x_limits(1) & x <= x_limits(2));
iy = find(y >= y_limits(1) & y <= y_limits(2));
assert(~isempty(ix) && ~isempty(iy), 'The requested crop is outside the data.');

xc = x(ix);
yc = y(iy);
e3c = e3(iy, ix);
nec = ne(iy, ix);

%% Identify tilted optical lobes from the two-dimensional zero-field pattern
% Estimate the dominant normal of the E3=0 contours from the intensity-
% weighted field-gradient tensor. Optical lobes are then separated by sign
% changes in the projected coordinate normal to those contours.
dx = mean(diff(xc));
dy = mean(diff(yc));
field_display = conv2(e3c, [0.2 0.6 0.2], 'same');
[dE_dx, dE_dy] = gradient(field_display, dx, dy);
weight = e3c.^2;
weight(weight < 0.03*max(weight(:))) = 0;
tensor = [sum(weight(:).*dE_dx(:).^2), sum(weight(:).*dE_dx(:).*dE_dy(:)); ...
          sum(weight(:).*dE_dx(:).*dE_dy(:)), sum(weight(:).*dE_dy(:).^2)];
[vectors, values] = eig(tensor);
[~, principal] = max(diag(values));
wave_normal = vectors(:, principal);
if wave_normal(1) < 0, wave_normal = -wave_normal; end
wave_tangent = [-wave_normal(2); wave_normal(1)];
disk_angle_xy = atan2(wave_normal(2), wave_normal(1));

[Xc, Yc] = meshgrid(xc, yc);
S = wave_normal(1)*(Xc-x_centroid) + wave_normal(2)*(Yc-y_centroid);
Q = wave_tangent(1)*(Xc-x_centroid) + wave_tangent(2)*(Yc-y_centroid);
q_center = sum(Q(:).*weight(:))/max(sum(weight(:)), eps);
projection_weight = exp(-0.5*((Q-q_center)/35).^2);

ds = 0.08;
s_axis = min(S(:)):ds:max(S(:));
s_bin = max(1, min(numel(s_axis), round((S-min(S(:)))/ds)+1));
projected_field = accumarray(s_bin(:), ...
    field_display(:).*projection_weight(:), [numel(s_axis),1], @sum, 0) ./ ...
    max(accumarray(s_bin(:), projection_weight(:), ...
    [numel(s_axis),1], @sum, 0), eps);
projected_field = conv(projected_field, [0.2 0.6 0.2], 'same');
signal_limit = 0.004*max(abs(projected_field));
projected_sign = sign(projected_field);
projected_sign(abs(projected_field)<signal_limit) = 0;
for k = 2:numel(projected_sign)
    if projected_sign(k)==0, projected_sign(k)=projected_sign(k-1); end
end
for k = numel(projected_sign)-1:-1:1
    if projected_sign(k)==0, projected_sign(k)=projected_sign(k+1); end
end

lobe_edges = [1; find(projected_sign(1:end-1).*projected_sign(2:end)<0)+1; ...
    numel(projected_sign)+1];
lobe_peak = zeros(numel(lobe_edges)-1,1);
lobe_energy = zeros(numel(lobe_edges)-1,1);
for k = 1:numel(lobe_peak)
    segment = lobe_edges(k):(lobe_edges(k+1)-1);
    [~, local_peak] = max(abs(projected_field(segment)));
    lobe_peak(k) = segment(local_peak);
    s_left_test = s_axis(lobe_edges(k));
    s_right_test = s_axis(lobe_edges(k+1)-1);
    lobe_energy(k) = sum(weight(S>=s_left_test & S<=s_right_test));
end
valid = abs(projected_field(lobe_peak)) >= signal_limit & ...
    lobe_energy >= 0.002*max(lobe_energy);
lobe_peak = lobe_peak(valid);
lobe_id_all = find(valid);
[~, order] = sort(s_axis(lobe_peak), 'descend');
lobe_peak = lobe_peak(order);
lobe_id_all = lobe_id_all(order);
assert(numel(lobe_peak)>=number_of_disks, ...
    'Only %d tilted optical lobes were found.', numel(lobe_peak));
lobe_peak = lobe_peak(1:number_of_disks);
lobe_id_all = lobe_id_all(1:number_of_disks);

disk_x = zeros(1,number_of_disks);
disk_y = zeros(1,number_of_disks);
disk_radius_y = zeros(1,number_of_disks);
disk_half_thickness = zeros(1,number_of_disks);
for k = 1:number_of_disks
    lobe_id = lobe_id_all(k);
    s_left = s_axis(lobe_edges(lobe_id));
    s_right = s_axis(lobe_edges(lobe_id+1)-1);
    lobe_mask = S>=s_left & S<=s_right;
    lobe_weight = weight.*lobe_mask;
    norm_weight = max(sum(lobe_weight(:)),eps);
    disk_x(k) = sum(Xc(:).*lobe_weight(:))/norm_weight;
    disk_y(k) = sum(Yc(:).*lobe_weight(:))/norm_weight;
    q_local = wave_tangent(1)*(Xc-disk_x(k)) + ...
        wave_tangent(2)*(Yc-disk_y(k));
    sigma_q = sqrt(sum(q_local(:).^2.*lobe_weight(:))/norm_weight);
    disk_radius_y(k) = min(max(1.8*sigma_q,3.0),18.0);
    disk_half_thickness(k) = min(maximum_disk_half_width, ...
        max(0.06,0.45*abs(s_right-s_left)));
end

%% Figure 1: 3D visualization of the original 2D PIC result
fig1 = figure(1);
set(fig1, 'Color', 'w', 'Position', [120 80 980 760]);
ax1 = axes(fig1);
hold(ax1, 'on');

% Downsample only for rendering; the physical crop and selected disk
% positions are calculated from the full-resolution arrays.
xs = xc(1:surface_stride:end);
ys = yc(1:surface_stride:end);
nes = nec(1:surface_stride:end, 1:surface_stride:end);

% Suppress grid-scale PIC noise only in the height-encoded rendering. The
% source arrays and all laser-position calculations remain unfiltered.
kernel_y = exp(-((-6:6)/3).^2/2);
kernel_x = exp(-((-18:18)/8).^2/2);
display_kernel = kernel_y(:) * kernel_x(:)';
display_kernel = display_kernel / sum(display_kernel(:));
nes = conv2(nes, display_kernel, 'same');

ne_scale = robust_upper_limit(nes, 0.995);
ne_normalized = min(nes ./ max(ne_scale, eps), 1);
[Xs, Ys] = meshgrid(xs, ys);

% Estimate the on-axis background density n0 upstream of the pulse. The
% height map is centered on this value, so ne=n0 and the laser centroid are
% both located at z=0. Height represents density variation, not a third
% simulated spatial coordinate.
reference_x = xs <= (x_centroid-20);
reference_y = abs(ys-y_centroid) <= 3;
reference_values = nes(reference_y, reference_x);
if isempty(reference_values)
    reference_values = nes(:);
end
n0 = median(reference_values(:));
ne_deviation = nes - n0;
deviation_scale = robust_upper_limit(abs(ne_deviation), 0.99);
Zs = ne_amplitude .* ne_deviation ./ max(deviation_scale, eps);
density_rgb = density_truecolor(ne_normalized);
% Remove the numerical-density floor geometrically. Relying on alpha alone
% leaves a tinted rectangular plate in MATLAB's exported raster image.
density_mask = ne_normalized >= 0.025;
Zs(~density_mask) = NaN;
density_alpha = ne_alpha .* min(ne_normalized./0.14, 1);
surf(ax1, Xs, Ys, Zs, density_rgb, ...
    'FaceColor', 'texturemap', 'EdgeColor', 'none', ...
    'AlphaData', density_alpha, 'AlphaDataMapping', 'none', ...
    'FaceAlpha', 'texturemap', 'DiffuseStrength', 0.72, ...
    'SpecularStrength', 0.08);

% Render the selected optical crests as finite-thickness disks. Their x
% positions, transverse centroids, and radii all come from the 2D e3 data.
blue = [0.04 0.28 0.95];
deep_blue = [0.01 0.05 0.28];
red = [0.92 0.06 0.03];
yellow = [1.00 0.62 0.02];
blue_drawn = 0;
red_drawn = 0;
for k = 1:numel(disk_x)
    if ~ismember(k, disk_visualized)
        continue;
    end
    if disk_x(k) > disk_split_um
        blue_drawn = blue_drawn + 1;
        if mod(blue_drawn,2) == 1, disk_color = blue; else, disk_color = deep_blue; end
    else
        red_drawn = red_drawn + 1;
        if mod(red_drawn,2) == 1, disk_color = red; else, disk_color = yellow; end
    end
    local_z_radius = max(laser_disk_z_radius, 0.48*disk_radius_y(k));
    local_x_half_width = disk_longitudinal_fatness * ...
        disk_half_thickness(k);
    draw_laser_disk(ax1, disk_x(k), disk_y(k), laser_plane_z, ...
        local_x_half_width, disk_radius_y(k), local_z_radius, ...
        disk_angle_xy, disk_color, 0.82);
end

view(ax1, fig1_view(1), fig1_view(2));
camup(ax1, [0 1 0]);
daspect(ax1, [1 1 0.72]);
camproj(ax1, 'perspective');
axis(ax1, 'tight');
axis(ax1, 'off');
set(ax1, 'FontName', 'Times New Roman', 'FontSize', 15, ...
    'Clipping', 'on');
lighting(ax1, 'gouraud');
camlight(ax1, 'headlight');


%% Figure 2: source 2D laser field and the six selected disk positions
if if_2d_visual
    fig2 = figure(2);
    set(fig2, 'Color', 'w', 'Position', [1150 120 920 520]);
    ax2 = axes(fig2);
    imagesc(ax2, xc, yc, e3c);
    set(ax2, 'YDir', 'normal');
    axis(ax2, 'tight');
    xlim(ax2, [min(disk_x)-2, max(disk_x)+2]);
    hold(ax2, 'on');
    cmax = robust_upper_limit(abs(e3c), 0.995);
    caxis(ax2, [-cmax cmax]);
    colormap(ax2, blue_white_red(256));

    for k = 1:numel(disk_x)
        if ~ismember(k, disk_visualized)
            continue;
        end
        line_coordinate = [-disk_radius_y(k), disk_radius_y(k)];
        plot(ax2, disk_x(k)+wave_tangent(1)*line_coordinate, ...
            disk_y(k)+wave_tangent(2)*line_coordinate, ...
            '--', 'Color', [0.05 0.05 0.05], 'LineWidth', 0.9);
        label_x = disk_x(k)-wave_tangent(1)*(disk_radius_y(k)+5);
        label_y = disk_y(k)-wave_tangent(2)*(disk_radius_y(k)+5);
        text(ax2, label_x, label_y, ...
            sprintf('%d', k), 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', 'FontName', ...
            'Times New Roman', 'FontSize', 10, 'FontWeight', 'bold');
    end

    xlabel(ax2, 'x (\mum)', 'FontName', 'Times New Roman');
    ylabel(ax2, 'y (\mum)', 'FontName', 'Times New Roman');
    set(ax2, 'FontName', 'Times New Roman', 'FontSize', 14, ...
        'LineWidth', 1.0, 'Layer', 'top');
    cb = colorbar(ax2);
    cb.Label.String = 'E_3/(m_ec\omega_0/e)';
    cb.Label.FontName = 'Times New Roman';
    title(ax2, '2D PIC field and selected thin-disk positions', ...
        'FontName', 'Times New Roman', 'FontWeight', 'normal');
end

%% Save editable and high-resolution outputs
savefig(fig1, fullfile(code_dir, 'Fig1_d_pseudo3D.fig'));
print(fig1, fullfile(code_dir, 'Fig1_d_pseudo3D.png'), '-dpng', '-r400');
if if_2d_visual
    savefig(fig2, fullfile(code_dir, 'Fig1_d_2D_visual.fig'));
    print(fig2, fullfile(code_dir, 'Fig1_d_2D_visual.png'), '-dpng', '-r400');
end

fprintf('Laser centroid: x = %.3f um, y = %.3f um\n', ...
    x_centroid, y_centroid);
fprintf('Selected disk positions (um):\n');
disp([disk_x(:), disk_y(:), disk_radius_y(:)]);

%% Local functions
function value = robust_upper_limit(data, fraction)
values = sort(data(isfinite(data)));
if isempty(values)
    value = 1;
    return;
end
index = max(1, min(numel(values), round(fraction*numel(values))));
value = values(index);
end

function rgb = density_truecolor(normalized_density)
% White-to-sky-blue density palette. Height remains the primary encoding.
% in wangzhan_fig.m. Height remains the primary density encoding.
low = [1.00 1.00 1.00];
mid = [0 191 255]/255;
high = [0 0 255]/255;
rgb = zeros([size(normalized_density), 3]);
for channel = 1:3
    first = low(channel) + ...
        (mid(channel)-low(channel)) .* min(2*normalized_density, 1);
    second_weight = max(2*normalized_density-1, 0);
    rgb(:,:,channel) = first + ...
        (high(channel)-mid(channel)) .* second_weight;
end
end

function draw_laser_disk(ax, x0, y0, z0, half_thickness, ...
        radius_y, radius_z, angle_xy, color, alpha_value)
% A longitudinally flattened ellipsoid avoids the hard cylindrical side
% wall while retaining the fitted width of each optical lobe.
[unit_x, unit_y, unit_z] = sphere(72);
x_local = half_thickness .* unit_x;
y_local = radius_y .* unit_y;
x_surface = x0 + cos(angle_xy).*x_local - sin(angle_xy).*y_local;
y_surface = y0 + sin(angle_xy).*x_local + cos(angle_xy).*y_local;
z_surface = z0 + radius_z .* unit_z;
surf(ax, x_surface, y_surface, z_surface, ...
    'FaceColor', color, 'EdgeColor', 'none', ...
    'FaceAlpha', alpha_value, 'DiffuseStrength', 0.70, ...
    'SpecularStrength', 0.25);
end

function cmap = blue_white_red(n)
if nargin < 1
    n = 256;
end
half = floor(n/2);
blue = [0.02 0.12 0.70];
white = [1 1 1];
red = [0.82 0.04 0.04];
cmap = [linspace(blue(1),white(1),half)', ...
        linspace(blue(2),white(2),half)', ...
        linspace(blue(3),white(3),half)'; ...
        linspace(white(1),red(1),n-half)', ...
        linspace(white(2),red(2),n-half)', ...
        linspace(white(3),red(3),n-half)'];
end
