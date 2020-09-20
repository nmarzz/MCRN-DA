% This script animates the height field and the vorticity produced by
% a shallow water model. It should be called only after shallow_water_model
% has been run.
clear all; close all; clc;
load('SWE_run_4days.mat');
modeloutput= x_save;
%% POD
% load('SWE_POD_r2.mat')
% load('SWE_POD_r50_new.mat')
load('SWE_POD_r30_new.mat')
modeloutput_POD=modeloutput_truncated;
%% DMD
% % load('SWE_DMD_r2.mat')
% load('SWE_DMD_r10.mat')
% % load('SWE_DMD_r20.mat')
% modeloutput=X_dmd;

L=length(t_save);%number of snapshots
m=length(y);%y-dimension
n=length(x);%x-dimension
N_gridpoints=m*n;%each snapshot has x*y-dimensions

u=modeloutput(1:N_gridpoints,:);v=modeloutput(N_gridpoints+1:N_gridpoints*2,:);h=modeloutput((N_gridpoints*2)+1:end,:);
u_save=reshape(u,n,m,L);v_save=reshape(v,n,m,L);h_save=reshape(h,n,m,L);

%%
u_POD=modeloutput_POD(1:N_gridpoints,:);
v_POD=modeloutput_POD(N_gridpoints+1:N_gridpoints*2,:);
h_POD=modeloutput_POD((N_gridpoints*2)+1:end,:);
u_save_POD=reshape(u_POD,n,m,L);
v_save_POD=reshape(v_POD,n,m,L);
h_save_POD=reshape(h_POD,n,m,L);
%%
h = squeeze(h_save(:,:,L));
u = squeeze(u_save(:,:,L));
v = squeeze(v_save(:,:,L));
%
h_POD = squeeze(h_save_POD(:,:,L));
u_POD = squeeze(u_save_POD(:,:,L));
v_POD= squeeze(v_save_POD(:,:,L));

interval=6;
x_1000km = x.*1e-6;
y_1000km = y.*1e-6;

colormap(cmocean('curl'));
slice=4000;
figure(9)
t = tiledlayout(2,1,'TileSpacing','Compact');
nexttile
%plot original velocity
h1=quiver(x_1000km(3:interval:end), y_1000km(3:interval:end), ...
    u(3:interval:end, 3:interval:end)',...
    v(3:interval:end, 3:interval:end)');
% h1=quiver( x, y, u(:,:,slice),v(:,:,slice));
set(h1,'AutoScale','on', 'AutoScaleFactor', 2,'color',[0 0 1])
set(gca,'xlim',[0 25],'Ylim',[0 5])
title('Truth')
nexttile
%plot truncated velocity
h2=quiver(x_1000km(3:interval:end), y_1000km(3:interval:end), ...
    u_POD(3:interval:end, 3:interval:end)',...
    v_POD(3:interval:end, 3:interval:end)');
set(h2,'AutoScale','on', 'AutoScaleFactor', 2,'color',[1 0 1])
set(gca,'xlim',[0 25],'Ylim',[0 5])
title('POD(r=50)')
%%
Error=modeloutput-modeloutput_POD;%difference matrix between X and X_r
u_E_reshaped=Error(1:N_gridpoints,:);
v_E_reshaped=Error(N_gridpoints+1:N_gridpoints*2,:);
h_E_reshaped=Error((N_gridpoints*2)+1:end,:);

Ux_E=reshape(u_E_reshaped,n,m,L);%Horizontal velocity of the difference matrix
Uy_E=reshape(v_E_reshaped,n,m,L);%Vertical velocity of the difference matrix
Uh_E=reshape(h_E_reshaped,n,m,L);%Vertical velocity of the difference matrix
u_EE = squeeze(Ux_E(:,:,L));
v_EE= squeeze(Uy_E(:,:,L));
h_EE= squeeze(Uy_E(:,:,L));
%%
figure(10)
TOLC=ptc12(9);
f=tiledlayout(2,1,'TileSpacing','Compact');
nexttile
%plot original velocity
h1=quiver(x_1000km(3:interval:end), y_1000km(3:interval:end), ...
    u(3:interval:end, 3:interval:end)',...
    v(3:interval:end, 3:interval:end)');

set(h1,'AutoScale','on', 'AutoScaleFactor', 2,'Color', TOLC(1,:))
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
hold on
%plot truncated velocity
h2=quiver(x_1000km(3:interval:end), y_1000km(3:interval:end), ...
    u_POD(3:interval:end, 3:interval:end)',...
    v_POD(3:interval:end, 3:interval:end)');
set(h2,'AutoScale','on', 'AutoScaleFactor', 2,'Color', TOLC(7,:))
set(gca,'xlim',[0 25],'Ylim',[0 5])
title('Comparison')
nexttile
%plot the error difference
hh=quiver(x_1000km(3:interval:end), y_1000km(3:interval:end), ...
    u_EE(3:interval:end, 3:interval:end)',...
    v_EE(3:interval:end, 3:interval:end)');
set(hh,'AutoScale','on', 'AutoScaleFactor', 2,'Color', TOLC(1,:))
set(gca,'xlim',[0 25],'Ylim',[0 5])
title('Difference')
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
%%
plot_height_range = [9500 10500];
if mean(plot_height_range) > 1000
    height_scale = 0.001;
    height_title = 'Height (km)';
else
    height_scale = 1;
    height_title = 'Height (m)';
end

for it = slice
    figure(11)
    t1 = tiledlayout(2,1,'TileSpacing','Compact');
    nexttile
    %plot original velocity
    handle = image(x_1000km, y_1000km, (h'+H').*height_scale);
    set(handle,'CDataMapping','scaled');
    hold on
    warning off
    contour(x_1000km, y_1000km, H',[1:1000:8001],'k');
    warning on
    
    quiver(x_1000km(3:interval:end), y_1000km(3:interval:end), ...
        u(3:interval:end, 3:interval:end)',...
        v(3:interval:end, 3:interval:end)','k');
    set(gca,'xlim',[0 25],'Ylim',[0 5])
    title('Model')
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
    nexttile
    %plot truncated velocity
    handle = image(x_1000km, y_1000km, (h_POD'+H').*height_scale);
    set(handle,'CDataMapping','scaled');
    hold on
    warning off
    contour(x_1000km, y_1000km, H',[1:1000:8001],'k');
    warning on
    
    quiver(x_1000km(3:interval:end), y_1000km(3:interval:end), ...
        u_POD(3:interval:end, 3:interval:end)',...
        v_POD(3:interval:end, 3:interval:end)','k');
    set(gca,'xlim',[0 25],'Ylim',[0 5])
    title('POD(r=30)')
    colormap(cmocean('curl'));
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
end
for it = slice
    figure(12)
    tiledlayout(2,1,'TileSpacing','Compact');
    nexttile;
    %plot original velocity
    handle = image(x_1000km, y_1000km, (h'+H').*height_scale);
    set(handle,'CDataMapping','scaled');
    hold on
    warning off
    contour(x_1000km, y_1000km, H',[1:1000:8001],'k');
    warning on
    
    quiver(x_1000km(3:interval:end), y_1000km(3:interval:end), ...
        u(3:interval:end, 3:interval:end)',...
        v(3:interval:end, 3:interval:end)','k');
    set(gca,'xlim',[0 25],'Ylim',[0 5])
    title('Model')
     set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
     caxis([-1 1])
    nexttile;
    %plot truncated velocity
    handle = image(x_1000km, y_1000km, (h_EE'+H').*height_scale);
    set(handle,'CDataMapping','scaled');
    hold on
    warning off
    contour(x_1000km, y_1000km, H',[1:1000:8001],'k');
    warning on
    
    quiver(x_1000km(3:interval:end), y_1000km(3:interval:end), ...
        u_EE(3:interval:end, 3:interval:end)',...
        v_EE(3:interval:end, 3:interval:end)','k');
    set(gca,'xlim',[0 25],'Ylim',[0 5])
    title('POD(r=30)')
    colormap('redblue')
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
    hold off
    colorbar('southoutside')



end