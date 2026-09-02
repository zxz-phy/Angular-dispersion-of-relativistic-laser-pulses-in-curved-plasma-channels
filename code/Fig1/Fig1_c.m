clc;
clear;
close all;

%% User controls
t = 3100;                    % OSIRIS output time used in the main text (um)
if_2d_visual = true;         % Figure 2: source 2D field with disk indices
number_of_disks = 18;
disk_visualized = [7:2:14 15:18];
disk_split_um = 3044.67;     % x>split: blue; x<split: red/yellow
leading_search_limit_um = 3055;
merged_lobe_center_um = 3042.77;
merged_lobe_range_um = [3041.0 3044.67];
trailing_search_limit_um = 3040.0; % reject the weak false lobe near 3041 um

% Display window around the laser centroid. Empty values trigger automatic
% limits around the intensity-weighted centroid.
x_window_um = 70;
y_window_um = 50;

% Pseudo-3D display parameters. The third coordinate is a visualization
% coordinate and must not be interpreted as a simulated spatial dimension.
laser_ne_offset = 5;         % laser center above the ne=n0 plane
ne_amplitude = 10;           % vertical gain of (ne-n0) in Figure 1
laser_plane_z = laser_ne_offset;
laser_disk_z_radius = 3.4;
maximum_disk_half_width = 0.75;
disk_longitudinal_fatness = 2; % multiplier along propagation direction x
ne_alpha = 0.90;

% Figure-1 box limits [min max] in um for x, y, and pseudo-z, and view
% angles [azimuth elevation] in degrees.
fig1_box_um = [3016 3064; ...
               -24   24; ...
               -12   20];
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

%% Select one thin disk from each signed optical lobe
% A transverse maximum of |E3| can contain several local maxima without a
% field reversal. Instead, construct a signed on-axis signal, divide it at
% genuine sign changes, and retain one representative crest per lobe.
axis_band = abs(yc-y_centroid) <= 1.0;
signed_signal = mean(e3c(axis_band, :), 1);
signed_signal = conv(signed_signal, [0.2 0.6 0.2], 'same');
signal_limit = 0.003 * max(abs(signed_signal));

field_sign = sign(signed_signal);
field_sign(abs(signed_signal) < signal_limit) = 0;
% Fill isolated near-zero samples without creating artificial sign flips.
for k = 2:numel(field_sign)
    if field_sign(k) == 0
        field_sign(k) = field_sign(k-1);
    end
end
for k = numel(field_sign)-1:-1:1
    if field_sign(k) == 0
        field_sign(k) = field_sign(k+1);
    end
end

lobe_edges = [1, find(field_sign(1:end-1).*field_sign(2:end) < 0)+1, ...
    numel(field_sign)+1];
candidate = zeros(1, 0);
for k = 1:numel(lobe_edges)-1
    segment = lobe_edges(k):(lobe_edges(k+1)-1);
    [segment_peak, local_index] = max(abs(signed_signal(segment)));
    if segment_peak >= signal_limit
        candidate(end+1) = segment(local_index); %#ok<SAGROW>
    end
end

candidate = candidate(xc(candidate) <= leading_search_limit_um);

% Treat the broad weak-field structure as one lobe and place its displayed
% center at the measured center rather than at an edge-local maximum.
in_merged_lobe = xc(candidate) >= merged_lobe_range_um(1) & ...
    xc(candidate) <= merged_lobe_range_um(2);
candidate(in_merged_lobe) = [];
[~, merged_center_index] = min(abs(xc-merged_lobe_center_um));
candidate = sort(candidate, 'descend');

leading_candidate = candidate(xc(candidate) > disk_split_um);
trailing_candidate = candidate(xc(candidate) < trailing_search_limit_um);
assert(numel(leading_candidate) >= 14 && numel(trailing_candidate) >= 3, ...
    'Insufficient signed lobes in the requested leading/trailing ranges.');
leading_pick = unique(round(linspace(1, numel(leading_candidate), 14)), ...
    'stable');
if numel(leading_pick) > 14
    leading_pick = leading_pick(1:14);
end
selected = [leading_candidate(leading_pick), merged_center_index, ...
    trailing_candidate(1:3)];

disk_x = xc(selected);
disk_y = zeros(size(selected));
disk_radius_y = zeros(size(selected));
disk_half_thickness = zeros(size(selected));
for k = 1:numel(selected)
    profile = e3c(:, selected(k)).^2;
    profile_sum = max(sum(profile), eps);
    disk_y(k) = sum(yc(:) .* profile) / profile_sum;
    sigma_y = sqrt(sum((yc(:)-disk_y(k)).^2 .* profile) / profile_sum);
    disk_radius_y(k) = min(max(1.8*sigma_y, 3.0), 14.0);

    % Fit the longitudinal thickness from the signed lobe containing this
    % crest. The rendered width therefore follows the local optical node.
    if abs(disk_x(k)-merged_lobe_center_um) <= mean(diff(xc))
        fitted_width = 0.45*diff(merged_lobe_range_um);
    else
        lobe_id = find(selected(k) >= lobe_edges(1:end-1) & ...
            selected(k) < lobe_edges(2:end), 1, 'first');
        left_index = lobe_edges(lobe_id);
        right_index = lobe_edges(lobe_id+1)-1;
        fitted_width = 0.45*abs(xc(right_index)-xc(left_index));
    end
    disk_half_thickness(k) = min(maximum_disk_half_width, ...
        max(0.06, fitted_width));
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
        disk_color, 0.82);
end

view(ax1, fig1_view(1), fig1_view(2));
camup(ax1, [0 1 0]);
axis(ax1, 'off');
daspect(ax1, [1 1 0.72]);
camproj(ax1, 'perspective');
set(ax1, 'XLim', fig1_box_um(1,:), 'YLim', fig1_box_um(2,:), ...
    'ZLim', fig1_box_um(3,:), 'XLimMode', 'manual', ...
    'YLimMode', 'manual', 'ZLimMode', 'manual');
camtarget(ax1, mean(fig1_box_um, 2)');
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
    ylim(ax2, [y_centroid-30, y_centroid+40]);
    hold(ax2, 'on');
    cmax = robust_upper_limit(abs(e3c), 0.995);
    caxis(ax2, [-cmax cmax]);
    colormap(ax2, blue_white_red(256));

    for k = 1:numel(disk_x)
        if ~ismember(k, disk_visualized)
            continue;
        end
        plot(ax2, [disk_x(k) disk_x(k)], ...
            [disk_y(k)-disk_radius_y(k), disk_y(k)+disk_radius_y(k)], ...
            '--', 'Color', [0.05 0.05 0.05], 'LineWidth', 0.9);
        label_y = y_centroid + 25 + 5*mod(k-1, 3);
        text(ax2, disk_x(k), label_y, ...
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
savefig(fig1, fullfile(code_dir, 'Fig1_c_pseudo3D.fig'));
print(fig1, fullfile(code_dir, 'Fig1_c_pseudo3D.png'), '-dpng', '-r400');
if if_2d_visual
    savefig(fig2, fullfile(code_dir, 'Fig1_c_2D_visual.fig'));
    print(fig2, fullfile(code_dir, 'Fig1_c_2D_visual.png'), '-dpng', '-r400');
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
        radius_y, radius_z, color, alpha_value)
% A longitudinally flattened ellipsoid avoids the hard cylindrical side
% wall while retaining the fitted width of each optical lobe.
[unit_x, unit_y, unit_z] = sphere(72);
x_surface = x0 + half_thickness .* unit_x;
y_surface = y0 + radius_y .* unit_y;
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
