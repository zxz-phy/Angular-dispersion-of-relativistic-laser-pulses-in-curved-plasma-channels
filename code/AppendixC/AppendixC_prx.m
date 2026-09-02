clc; clear; close all;

if_export = true;
fs = 20;

script_dir = fileparts(mfilename('fullpath'));
figs_root = fileparts(fileparts(script_dir));
data_dir = fullfile(figs_root, 'data', 'AppendixC');
output_dir = fullfile(figs_root, 'Figs', 'AppendixC');
if if_export && ~isfolder(output_dir)
    mkdir(output_dir);
end
%%
fg=figure(22);
fg.Position=[395.67       89.667       944.67       758.67];
%%
load(fullfile(data_dir, 'discharge_spectrum_20psi.mat'), ...
    'tdata', 'mymap', 'Ax', 'capx', 'A', 'f1', 'ne0', 'r0');
s1=subplot(2,2,1);
imagesc([587.6-19,587.6+19],[-120*1.6,120*1.6],tdata/max(tdata(:)))
colormap(mymap)
xlabel('\lambda (nm)')
ylabel('r (\mum)')
cb=colorbar;
cb.Location='eastoutside';
cb.Label.String='Intensity (a. u.)';
cb.Label.Rotation=90;
cb.Label.Clipping='off';
cb.FontSize=fs-2;
cb.Label.FontSize=fs-2;
caxis([0 1])
axis([587.6-10 587.6+10 -200 200])
hold on
title(['20 Psi'])
set(gca,'Ydir','Normal')
add_panel_label(s1, '(a)', fs)
s2=subplot(2,2,2);
hold on
dd=5;
for j=1:ceil(length(Ax)/dd)
    Ax2(j)=mean(Ax(dd*(j-1)+1:min(dd*j,length(Ax))));
end
cxx2=linspace(-120*1.6,120*1.6,length(Ax2));
% plot(cxx2,Ax2,'b.')
plot(capx,A,'b.')
p2=plot(f1,'r-');
p2.LineWidth=1;
lgd2=legend({'Plasma density n_e(r)', ...
    sprintf('n_0=%.2f \\times10^{18} cm^{-3}; r_0=%d \\mum', ...
    ne0/1e18, round(r0))}, 'Interpreter', 'tex')
lgd2.FontSize=fs-4;
xlabel('r (\mum)')
ylabel('n_e (10^{18} cm^{-3})')
axis([-200 200 0 3])
add_panel_label(s2, '(b)', fs)

%%
load(fullfile(data_dir, 'discharge_spectrum_10psi.mat'), ...
    'tdata', 'mymap', 'Ax', 'cx', 'A', 'f1', 'ne0', 'r0');
s3=subplot(2,2,3);
imagesc([587.6-19,587.6+19],[-120*1.6,120*1.6],tdata/max(tdata(:))*2)
colormap(mymap)
axis([587.6-10 587.6+10 -200 200])
xlabel('\lambda (nm)')
ylabel('r (\mum)')
caxis([0 1])
hold on
title(['10 Psi'])
set(gca,'Ydir','Normal')
add_panel_label(s3, '(c)', fs)
s4=subplot(2,2,4);
hold on
dd=5;
for j=1:ceil(length(Ax)/dd)
    Ax2(j)=mean(Ax(dd*(j-1)+1:min(dd*j,length(Ax))));
end
% cxx2=linspace(-120*1.6,120*1.6,length(Ax2));
% cxx2=linspace(-150,150,length(Ax2));
plot(cx,A,'b.')
p2=plot(f1,'r-');
p2.LineWidth=1;
lgd4=legend({'Plasma density n_e(r)', ...
    sprintf('n_0=%.2f \\times10^{18} cm^{-3}; r_0=%d \\mum', ...
    ne0/1e18, round(r0))}, 'Interpreter', 'tex')
lgd4.FontSize=fs-2;
xlabel('r (\mum)')
ylabel('n_e (10^{18} cm^{-3})')
axis([-200 200 0 3])
add_panel_label(s4, '(d)', fs)


%%
cb.Position=[0.910      0.58348     0.018581      0.34116];
s1.Position=[0.15      0.58384      0.33466      0.34116];
s3.Position=[0.15         0.11      0.33466      0.34116];
apply_axes_font([s1, s2, s3, s4], fs)
lgd2.FontSize=fs-6;
lgd4.FontSize=fs-6;

if if_export
    exportgraphics(fg, fullfile(output_dir, 'FigC1.png'), 'Resolution', 600);
    exportgraphics(fg, fullfile(output_dir, 'FigC1.pdf'), 'ContentType', 'vector');
    print(fg, fullfile(output_dir, 'FigC1.eps'), '-depsc', '-painters');
    fprintf('Saved Fig. C1 to %s\n', output_dir);
else
    fprintf('Fig. C1 export skipped (if_export = false).\n');
end

function add_panel_label(ax, label_text, fs)
text(ax, -0.12, 1.05, label_text, 'Units', 'normalized', ...
    'FontWeight', 'bold', 'FontSize', fs+4, 'VerticalAlignment', 'bottom', ...
    'Color', 'k', 'Clipping', 'off');
end

function apply_axes_font(ax_array, fs)
for ax = ax_array
    ax.FontSize = fs;
    ax.XLabel.FontSize = fs;
    ax.YLabel.FontSize = fs;
    ax.Title.FontSize = fs;
end
end





