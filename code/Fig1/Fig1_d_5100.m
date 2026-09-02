clc;
clear;
close all;

%% User controls
t = 5100;                    % OSIRIS output time used in the main text (um)
if_2d_visual = true;         % Figure 2: source 2D field with disk indices
number_of_nir_disks = 20;
number_of_mir_disks = 10;
number_of_disks = number_of_nir_disks + number_of_mir_disks;
N_visualize = 2:10;         % NIR family-local disk indices to render
M_visualize = 3:8;          % MIR family-local disk indices to render
assert(all(ismember(N_visualize, 1:number_of_nir_disks)), ...
    'N_visualize contains an index outside the detected NIR family.');
assert(all(ismember(M_visualize, 1:number_of_mir_disks)), ...
    'M_visualize contains an index outside the detected MIR family.');
disk_visualized = [N_visualize, number_of_nir_disks + M_visualize];
% Select the coherent NIR packet and exclude its rapidly oscillating front.
nir_roi_um = [5008 5031 -68 -28];  % [xmin xmax ymin ymax]
% Isolate the coherent MIR packet rather than its weak downstream tail.
mir_roi_um = [4965 4997 -112 -78];
analysis_roi_um = [4950 5059 -155 -10];
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

% Figure-1 box limits [min max] in um for x, y, and pseudo-z, and view
% angles [azimuth elevation] in degrees.
fig1_box_um = [4950 5050; ...
               -120  -20; ...
                -12   20];
fig1_view = [0 90];
surface_stride = 4;         % render every fourth density sample in x and y

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

x_limits = analysis_roi_um(1:2);
y_limits = analysis_roi_um(3:4);
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

%% Identify NIR and MIR lobes independently inside two spatial ROIs
% Each pulse family obtains its own zero-field contour orientation. This is
% required because the NIR and MIR components follow different trajectories.
nir_lobes = identify_tilted_lobes(e3c, xc, yc, nir_roi_um, ...
    number_of_nir_disks, maximum_disk_half_width);
mir_lobes = identify_tilted_lobes(e3c, xc, yc, mir_roi_um, ...
    number_of_mir_disks, maximum_disk_half_width);

disk_x = [nir_lobes.x, mir_lobes.x];
disk_y = [nir_lobes.y, mir_lobes.y];
disk_radius_y = [nir_lobes.radius, mir_lobes.radius];
disk_half_thickness = [nir_lobes.half_width, mir_lobes.half_width];
disk_angle_xy = [nir_lobes.angle, mir_lobes.angle];
disk_tangent_x = [nir_lobes.tangent_x, mir_lobes.tangent_x];
disk_tangent_y = [nir_lobes.tangent_y, mir_lobes.tangent_y];
disk_group = [ones(1,number_of_nir_disks), ...
    2*ones(1,number_of_mir_disks)];
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
    if disk_group(k) == 1
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
        disk_angle_xy(k), disk_color, 0.82);
end

view(ax1, fig1_view(1), fig1_view(2));
camup(ax1, [0 1 0]);
daspect(ax1, [1 1 0.72]);
camproj(ax1, 'perspective');
axis(ax1, 'off');
set(ax1, 'XLim', fig1_box_um(1,:), 'YLim', fig1_box_um(2,:), ...
    'ZLim', fig1_box_um(3,:), 'XLimMode', 'manual', ...
    'YLimMode', 'manual', 'ZLimMode', 'manual');
camtarget(ax1, mean(fig1_box_um, 2)');
set(ax1, 'FontName', 'Times New Roman', 'FontSize', 15, ...
    'Clipping', 'on');
lighting(ax1, 'gouraud');
camlight(ax1, 'headlight');


%% Figure 2: source 2D laser field and the selected NIR/MIR lobes
if if_2d_visual
    fig2 = figure(2);
    set(fig2, 'Color', 'w', 'Position', [1150 120 920 520]);
    ax2 = axes(fig2);
    imagesc(ax2, xc, yc, e3c);
    set(ax2, 'YDir', 'normal');
    axis(ax2, 'tight');
    xlim(ax2, analysis_roi_um(1:2));
    ylim(ax2, analysis_roi_um(3:4));
    hold(ax2, 'on');
    cmax = robust_upper_limit(abs(e3c), 0.995);
    caxis(ax2, [-cmax cmax]);
    colormap(ax2, blue_white_red(256));

    for k = 1:numel(disk_x)
        if ~ismember(k, disk_visualized)
            continue;
        end
        % Keep the marker within the local pulse envelope.
        line_half_length = 0.72*disk_radius_y(k);
        line_coordinate = [-line_half_length, line_half_length];
        plot(ax2, disk_x(k)+disk_tangent_x(k)*line_coordinate, ...
            disk_y(k)+disk_tangent_y(k)*line_coordinate, ...
            '--', 'Color', [0.05 0.05 0.05], 'LineWidth', 0.9);

        % Put compact family-local labels beside the same end of each lobe.
        label_x = disk_x(k)-disk_tangent_x(k)*(line_half_length+1.2);
        label_y = disk_y(k)-disk_tangent_y(k)*(line_half_length+1.2);
        if disk_group(k) == 1
            label_string = sprintf('N%d', k);
        else
            label_string = sprintf('M%d', k-number_of_nir_disks);
        end
        text(ax2, label_x, label_y, ...
            label_string, 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', 'FontName', ...
            'Times New Roman', 'FontSize', 9, 'FontWeight', 'bold', ...
            'Clipping', 'on');
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
savefig(fig1, fullfile(code_dir, 'Fig1_d_5100_pseudo3D.fig'));
print(fig1, fullfile(code_dir, 'Fig1_d_5100_pseudo3D.png'), '-dpng', '-r400');
if if_2d_visual
    savefig(fig2, fullfile(code_dir, 'Fig1_d_5100_2D_visual.fig'));
    print(fig2, fullfile(code_dir, 'Fig1_d_5100_2D_visual.png'), '-dpng', '-r400');
end

fprintf('Laser centroid: x = %.3f um, y = %.3f um\n', ...
    x_centroid, y_centroid);
fprintf('Selected disk positions (um):\n');
disp([disk_x(:), disk_y(:), disk_radius_y(:)]);

%% Local functions
function lobes = identify_tilted_lobes(field, x, y, roi, requested_count, max_half_width)
ix = find(x>=roi(1) & x<=roi(2));
iy = find(y>=roi(3) & y<=roi(4));
assert(~isempty(ix) && ~isempty(iy), 'Empty pulse-family ROI.');
E = field(iy,ix);
xr = x(ix);
yr = y(iy);
[X,Y] = meshgrid(xr,yr);
dx = mean(diff(xr));
dy = mean(diff(yr));
E = conv2(E,[0.2 0.6 0.2],'same');
I = E.^2;
I(I<0.02*max(I(:))) = 0;
[gx,gy] = gradient(E,dx,dy);
J = [sum(I(:).*gx(:).^2),sum(I(:).*gx(:).*gy(:)); ...
     sum(I(:).*gx(:).*gy(:)),sum(I(:).*gy(:).^2)];
[V,D] = eig(J);
[~,id] = max(diag(D));
n = V(:,id);
if n(1)<0, n=-n; end
tangent = [-n(2);n(1)];
normI = max(sum(I(:)),eps);
x0 = sum(X(:).*I(:))/normI;
y0 = sum(Y(:).*I(:))/normI;
S = n(1)*(X-x0)+n(2)*(Y-y0);
Q = tangent(1)*(X-x0)+tangent(2)*(Y-y0);
q0 = sum(Q(:).*I(:))/normI;
W = exp(-0.5*((Q-q0)/30).^2);
ds = 0.06;
saxis = min(S(:)):ds:max(S(:));
ibin = max(1,min(numel(saxis),round((S-min(S(:)))/ds)+1));
den = accumarray(ibin(:),W(:),[numel(saxis),1],@sum,0);
proj = accumarray(ibin(:),E(:).*W(:),[numel(saxis),1],@sum,0)./max(den,eps);
proj = conv(proj,[0.2 0.6 0.2],'same');
limit = 0.004*max(abs(proj));
sgn = sign(proj);
sgn(abs(proj)<limit)=0;
for j=2:numel(sgn),if sgn(j)==0,sgn(j)=sgn(j-1);end,end
for j=numel(sgn)-1:-1:1,if sgn(j)==0,sgn(j)=sgn(j+1);end,end
edges = [1;find(sgn(1:end-1).*sgn(2:end)<0)+1;numel(sgn)+1];
peaks = zeros(numel(edges)-1,1);
energy = zeros(size(peaks));
for j=1:numel(peaks)
    seg=edges(j):(edges(j+1)-1);
    [~,p]=max(abs(proj(seg)));
    peaks(j)=seg(p);
    mask=S>=saxis(edges(j)) & S<=saxis(edges(j+1)-1);
    energy(j)=sum(I(mask));
end
valid=abs(proj(peaks))>=limit & energy>=0.002*max(energy);
ids=find(valid);
[~,ord]=sort(saxis(peaks(valid)),'descend');
ids=ids(ord);
assert(numel(ids)>=requested_count,'Only %d lobes found in ROI.',numel(ids));
ids=ids(1:requested_count);

lobes.x=zeros(1,requested_count);
lobes.y=zeros(1,requested_count);
lobes.radius=zeros(1,requested_count);
lobes.half_width=zeros(1,requested_count);
for j=1:requested_count
    left=saxis(edges(ids(j)));
    right=saxis(edges(ids(j)+1)-1);
    mask=S>=left & S<=right;
    WI=I.*mask;
    nw=max(sum(WI(:)),eps);
    lobes.x(j)=sum(X(:).*WI(:))/nw;
    lobes.y(j)=sum(Y(:).*WI(:))/nw;
    qlocal=tangent(1)*(X-lobes.x(j))+tangent(2)*(Y-lobes.y(j));
    sigma=sqrt(sum(qlocal(:).^2.*WI(:))/nw);
    lobes.radius(j)=min(max(1.8*sigma,3),20);
    lobes.half_width(j)=min(max_half_width,max(0.06,0.45*abs(right-left)));
end
lobes.angle=repmat(atan2(n(2),n(1)),1,requested_count);
lobes.tangent_x=repmat(tangent(1),1,requested_count);
lobes.tangent_y=repmat(tangent(2),1,requested_count);
end

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
