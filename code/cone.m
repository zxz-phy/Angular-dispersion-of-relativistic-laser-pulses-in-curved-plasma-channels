clc;
clear;
close all;

%% Cone controls
cone_length = 20;       % Cone height/axial length
cone_rad = 1/40*pi;     % Cone half-opening angle (rad)

%% Geometry: apex at x = 0, circular base at x = cone_length
n_axial = 160;
n_azimuth = 240;
x_axis = linspace(0, cone_length, n_axial);
azimuth = linspace(0, 2*pi, n_azimuth);
[x_surface, azimuth_surface] = meshgrid(x_axis, azimuth);
radius_surface = x_surface*tan(cone_rad);
y_surface = radius_surface.*cos(azimuth_surface);
z_surface = radius_surface.*sin(azimuth_surface);
base_radius = cone_length*tan(cone_rad);

%% Figure
fg = figure(1);
set(fg, 'Color', 'w', 'Units', 'pixels', 'Position', [280, 150, 900, 600]);
ax = axes(fg);
hold(ax, 'on');

cone_color = [0.88, 0.035, 0.025];
surf(ax, x_surface, y_surface, z_surface, ...
    'FaceColor', cone_color, ...
    'EdgeColor', 'none', ...
    'FaceLighting', 'gouraud', ...
    'AmbientStrength', 0.22, ...
    'DiffuseStrength', 0.72, ...
    'SpecularStrength', 0.95, ...
    'SpecularExponent', 35, ...
    'SpecularColorReflectance', 0.9);

% Close the circular base of the cone.
base_azimuth = linspace(0, 2*pi, n_azimuth);
patch(ax, cone_length*ones(size(base_azimuth)), ...
    base_radius*cos(base_azimuth), ...
    base_radius*sin(base_azimuth), cone_color, ...
    'EdgeColor', 'none', ...
    'FaceLighting', 'gouraud', ...
    'AmbientStrength', 0.22, ...
    'DiffuseStrength', 0.72, ...
    'SpecularStrength', 0.8, ...
    'SpecularExponent', 28);

%% Glossy lighting and top view
material(ax, 'shiny');
lighting(ax, 'gouraud');
camlight(ax, 'headlight');
light(ax, 'Position', [0.7*cone_length, -2*base_radius, 2.2*base_radius], ...
    'Style', 'local', 'Color', [1, 0.88, 0.82]);

axis(ax, 'equal');
axis(ax, 'off');
xlim(ax, [-0.08*cone_length, 1.08*cone_length]);
ylim(ax, 1.2*[-base_radius, base_radius]);
zlim(ax, 1.2*[-base_radius, base_radius]);

% The cone axis remains horizontal along x; look down from +z onto the x-y plane.
view(ax, [0, 90]);
camproj(ax, 'perspective');
set(ax, 'Clipping', 'off');
