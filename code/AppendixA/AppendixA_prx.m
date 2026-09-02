clc; clear; close all;

script_dir = fileparts(mfilename('fullpath'));
figs_root = fileparts(fileparts(script_dir));
data_dir = fullfile(figs_root, 'data', 'AppendixA');
output_dir = fullfile(figs_root, 'Figs', 'AppendixA');
if ~isfolder(output_dir), mkdir(output_dir); end

trajectory = load(fullfile(data_dir, 'total_data_cen_angle_recal.mat'));
profile = load(fullfile(data_dir, 'deltay.mat'));
spectrum = load(fullfile(data_dir, 'theta_lamda_L3R20_5200.mat'));
field_file = fullfile(data_dir, 'e3-005200.h5');
e3 = double(h5read(field_file, '/e3'))' / 7.854;
xmax = h5readatt(field_file, '/', 'XMAX'); xmin = h5readatt(field_file, '/', 'XMIN');
x_field = linspace(xmin(1), xmax(1), size(e3,2)) / 1e3;
y_field = linspace(xmin(2), xmax(2), size(e3,1));
% Crop to the displayed exit region before rendering. The source HDF5 remains intact.
x_keep = x_field >= 5.05 & x_field <= 5.15;
y_keep = y_field >= -150 & y_field <= 0;
e3_display = e3(y_keep, x_keep);
x_field = x_field(x_keep); y_field = y_field(y_keep);
clear e3 x_keep y_keep;

fs = 16; fg = figure(6);
set(fg, 'Color', 'w', 'Units', 'pixels', 'Position', [100, 80, 1050, 780]);

%% (a) Propagation trajectories and sampled channel density
ax_a = axes(fg, 'Position', [0.10, 0.61, 0.82, 0.31]); hold(ax_a, 'on');
straight_length=3; curved_length=2; radius=20; channel_radius=.2;
x_curve=linspace(0,curved_length,2000); y_curve=sqrt(radius^2-x_curve.^2)-radius;
x_axis=[linspace(0,straight_length,3000), x_curve+straight_length];
y_axis=[zeros(1,3000), y_curve];
h_wall = plot(ax_a,x_axis+.1,y_axis+channel_radius,'k-');
plot(ax_a,x_axis+.1,y_axis-channel_radius,'k-');
h_drive = plot(ax_a,trajectory.cen_else(:,2)/1e3,trajectory.cen_else(:,1)/1e3,'b-','LineWidth',2.5);
h_mir = plot(ax_a,trajectory.cen_ir(:,2)/1e3,trajectory.cen_ir(:,1)/1e3,'r-','LineWidth',2.5);
xlim(ax_a,[0 9.8]); ylim(ax_a,[-.8 .3]); box(ax_a,'on');
xlabel(ax_a,'$x~(\mathrm{mm})$','Interpreter','latex');
ylabel(ax_a,'$y~(\mathrm{mm})$','Interpreter','latex');
legend(ax_a,[h_wall,h_drive,h_mir],{'Channel wall','Drive pulse','MIR pulse'}, ...
    'Location','northeast','Box','off');
panel(ax_a,'(a)',fs);

ax_density=axes(fg,'Position',[0.102,0.615,0.43,0.095]);
ne=1+22.63*profile.deltay.^2/16^4; ne([1,end])=0;
plot(ax_density,profile.cenxs/1e3,ne,'k-','LineWidth',1.8); xlim(ax_density,[0 5.15]);
ylim(ax_density,[0 3]); set(ax_density,'XTick',[],'YTick',[1 2 3], ...
    'YAxisLocation','right','FontName','Times New Roman','FontSize',fs-2,'Box','off');
ylabel(ax_density,'$n_{\mathrm{wit}}/n_0$','Interpreter','latex');

%% (b) Exit field with channel-density envelope
ax_b = axes(fg, 'Position', [0.10, 0.12, 0.35, 0.38]);
imagesc(ax_b,x_field,y_field,e3_display); set(ax_b,'YDir','normal');
xlim(ax_b,[5.05 5.15]); ylim(ax_b,[-150 0]); caxis(ax_b,[-1 1]);
colormap(ax_b, blue_white_red(256));
xlabel(ax_b,'$x~(\mathrm{mm})$','Interpreter','latex');
ylabel(ax_b,'$y~(\mu\mathrm{m})$','Interpreter','latex'); panel(ax_b,'(b)',fs);
cb_b=colorbar(ax_b); cb_b.Label.String='a_0';

%% (c) Output angular spectrum
ax_c = axes(fg, 'Position', [0.57, 0.12, 0.35, 0.38]);
axes(ax_c); %#ok<LAXES>
% The stored coordinate grid is the publication-resolution 10x display grid.
% Sample the full FFT intensity onto that grid without altering the source MAT file.
Y_display = spectrum.Y(1:10:end, 1:10:end);
warp(spectrum.rho,spectrum.theta,spectrum.Z,Y_display); view(ax_c,2);
clear Y_display spectrum;
set(ax_c,'YDir','normal');
xlim(ax_c,[0 10]); ylim(ax_c,[-20 20]); caxis(ax_c,[0 8e3]);
colormap(ax_c, faded_jet(256));
xlabel(ax_c,'$\lambda~(\mu\mathrm{m})$','Interpreter','latex');
ylabel(ax_c,'$\theta~(\mathrm{deg})$','Interpreter','latex'); panel(ax_c,'(c)',fs);
cb_c=colorbar(ax_c); cb_c.Label.String='Intensity (a. u.)';

set([ax_a,ax_b,ax_c],'FontName','Times New Roman','FontSize',fs,'Box','on','TickDir','out');
exportgraphics(fg,fullfile(output_dir,'FigA1.png'),'Resolution',600);
print(fg,fullfile(output_dir,'FigA1.eps'),'-depsc','-painters');
exportgraphics(fg,fullfile(output_dir,'FigA1.pdf'),'ContentType','vector');
fprintf('Saved Fig. A1 to %s\n',output_dir);

function panel(ax,label,fs)
text(ax,.01,1.02,label,'Units','normalized','FontName','Times New Roman', ...
    'FontSize',fs+6,'FontWeight','bold','VerticalAlignment','bottom','Clipping','off');
end
function map=blue_white_red(n)
h=floor(n/2); map=[linspace(0,1,h)',linspace(0,1,h)',ones(h,1); ...
    ones(n-h,1),linspace(1,0,n-h)',linspace(1,0,n-h)'];
end
function map=faded_jet(n)
map=jet(n); cut=round(.2*n); map(1:cut,:)=[linspace(1,map(cut,1),cut)', ...
    linspace(1,map(cut,2),cut)',linspace(1,map(cut,3),cut)'];
end
