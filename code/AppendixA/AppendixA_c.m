clc; clear; close all;

% Resolve all dependencies relative to this script so the Figs folder is portable.
script_dir = fileparts(mfilename('fullpath'));
figs_dir = fileparts(fileparts(script_dir));
data_dir = fullfile(figs_dir, 'data', 'AppendixA');
output_dir = fullfile(figs_dir, 'Figs', 'AppendixA');

load(fullfile(data_dir, 'theta_lamda_L3R20_5200.mat'))
% load('theta_lamda_L3R20_3100_small_used.mat')
% axnum=10;
% load('theta_lamda_L3R20_5200_with_aperture_0429.mat')
axnum=20/180*pi;
%%
% blue1=[100 100 255]/255;
% blue2=[255 255 0]/255;
% mymap=[  linspace(blue1(1),blue2(1),1e3)',...
%          linspace(blue1(2),blue2(2),1e3)', ... 
%          linspace(blue1(3),blue2(3),1e3)'];
mymap=jet(512);
% mymap=mymap(1:420,:);
sj=size(mymap);
duct=round(sj(1)*0.25);
% duct=1; 
ms(:,:)=[linspace(1,mymap(duct,1),duct)',...
         linspace(1,(mymap(duct,2)),duct)', ... 
         linspace(1,(mymap(duct,3)),duct)'];
ms=ms.^(1/8); 
mymap=[ms;mymap(round((1+duct)):end,:)];
figure('Visible', 'off')
hold on
Y=Y/1e4;
% Convert the laboratory-frame deflection to the relative deviation angle.
% The sign inversion also restores the convention used in the manuscript.
theta0_deg = rad2deg(2/20);
theta = -theta - theta0_deg;
warp(rho,theta,Z,Y), view(2), axis square tight off
colormap(mymap)
set(gca,'Ydir','Normal')
caxis([0 1])
box on
%%
hold on
% plot([0 0],[-pi/6 pi/6],'k-','LineWidth',0.1)
% plot([0 15],[pi/6 pi/6],'k-','LineWidth',0.1)
% plot([15 15],[-pi/6 pi/6],'k-','LineWidth',0.1)
% plot([0 15],[-pi/6 -pi/6],'k-','LineWidth',0.1)
% offaxis_deg=5.41;
axis([0 10 -20 20])
axis on
grid on
box on
xup=6;
%%
% xdown=-7;
% plot([0 15],[xdown xdown],'k-','LineWidth',0.1)
% xdown=-7;
% plot([0 15],[xdown xdown],'k-','LineWidth',0.1)
xl=-0;
plot([xl xl],[-20 20],'k-','LineWidth',0.1)
xr=10;
plot([xr xr],[-20 20],'k-','LineWidth',0.1)
%%
% plot([0 0],[xdown xup],'k-','LineWidth',0.1)
% plot([0 15],[xup xup],'k-','LineWidth',0.1)
% plot([15 15],[xdown xup],'k-','LineWidth',0.1)

% axis([0 15  xdown xup])
axis off
grid off
box off

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
% OpenGL rasterizes the dense background inside the EPS and avoids an
% impractically large vector file from the high-resolution angular spectrum.
print(gcf, fullfile(output_dir, 'FigA1_C.eps'), '-depsc2', '-opengl', '-r600');
