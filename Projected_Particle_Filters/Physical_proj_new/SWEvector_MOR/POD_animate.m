% This script animates the height field and the vorticity produced by
% a shallow water model. It should be called only after shallow_water_model
% has been run.
clear all; close all; clc;
% load('SWE_run_1day.mat')
% modeloutput= x_save;
%% POD
% load('SWE_POD_r2.mat')
% load('SWE_POD_r10.mat')
load('SWE_POD_r20_new.mat')
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
slice=L;
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
title('POD(r=20)')
%%
Error=modeloutput-modeloutput_POD;%difference matrix between X and X_r
u_E_reshaped=Error(1:N_gridpoints,:);
v_E_reshaped=Error(N_gridpoints+1:N_gridpoints*2,:);
h_E_reshaped=Error((N_gridpoints*2)+1:end,:);

Ux_E=reshape(u_E_reshaped,n,m,L);%Horizontal velocity of the difference matrix
Uy_E=reshape(v_E_reshaped,n,m,L);%Vertical velocity of the difference matrix

u_EE = squeeze(Ux_E(:,:,L));
v_EE= squeeze(Uy_E(:,:,L));

%%
figure(10)
f=tiledlayout(2,1,'TileSpacing','Compact');
nexttile
%plot original velocity
h1=quiver(x_1000km(3:interval:end), y_1000km(3:interval:end), ...
    u(3:interval:end, 3:interval:end)',...
    v(3:interval:end, 3:interval:end)');

set(h1,'AutoScale','on', 'AutoScaleFactor', 2.5,'color',[0 0 1])
hold on
%plot truncated velocity
h2=quiver(x_1000km(3:interval:end), y_1000km(3:interval:end), ...
    u_POD(3:interval:end, 3:interval:end)',...
    v_POD(3:interval:end, 3:interval:end)');
set(h2,'AutoScale','on', 'AutoScaleFactor', 2,'color',[1 0 1])
set(gca,'xlim',[0 25],'Ylim',[0 5])
title('Comparison')
nexttile
%plot the error difference
hh=quiver(x_1000km(3:interval:end), y_1000km(3:interval:end), ...
    u_EE(3:interval:end, 3:interval:end)',...
    v_EE(3:interval:end, 3:interval:end)');
set(hh,'AutoScale','on', 'AutoScaleFactor', 2,'color',[0 0 1])
set(gca,'xlim',[0 25],'Ylim',[0 5])
title('Difference')
