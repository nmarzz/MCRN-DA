clear all; close all; clc;
load('SWE_run_4days.mat');
% modeloutput= x_save;
load('SWE_POD.mat')
modeloutput= diff_plot;
L=length(modeloutput(1,:));%number of snapshots
m=length(y);%y-dimension
n=length(x);%x-dimension
N_gridpoints=m*n;%each snapshot has x*y-dimensions

u=modeloutput(1:N_gridpoints,:);v=modeloutput(N_gridpoints+1:N_gridpoints*2,:);h=modeloutput((N_gridpoints*2)+1:end,:);
u_save=reshape(u,n,m,L);v_save=reshape(v,n,m,L);h_save=reshape(h,n,m,L);
x_1000km = x.*1e-6;
y_1000km = y.*1e-6;
plot_height_range = [9500 10500];
% Set colormap to have 64 entries
ncol=64;
% colormap(jet(ncol));
colormap(cmocean('curl'));
% Interval between arrows in the velocity vector plot
interval = 6;

% Set this to "true" to save each frame as a png file
plot_frames = false;

% Decide whether to show height in metres or km
if mean(plot_height_range) > 1000
    height_scale = 0.001;
    height_title = 'Velocity Field ($u$ and $v$), Difference';
else
    height_scale = 1;
    height_title = 'Velocity Field ($u$ and $v$), Difference';
end
% it=60*24;
it=50;
% Extract the height and velocity components for this frame
h = squeeze(h_save(:,:,it));
u = squeeze(u_save(:,:,it));
v = squeeze(v_save(:,:,it));
% Plot the height field
% it_1=60*48;
 it_1=100;
% Extract the height and velocity components for this frame
h_1= squeeze(h_save(:,:,it_1));
u_1 = squeeze(u_save(:,:,it_1));
v_1 = squeeze(v_save(:,:,it_1));
% it_2=60*72;
% % Extract the height and velocity components for this frame
% h_2= squeeze(h_save(:,:,it_2));
% u_2 = squeeze(u_save(:,:,it_2));
% v_2 = squeeze(v_save(:,:,it_2));
%   warning on
figure(1)
tt = tiledlayout(2,1);
hh(1)=nexttile;
qp=quiver(x_1000km(3:interval:end), y_1000km(3:interval:end), ...
    u(3:interval:end, 3:interval:end)',...
    v(3:interval:end, 3:interval:end)');
set(qp,'AutoScale','on', 'AutoScaleFactor', 2,'color',[0 0 0])
legend('Scale: 2m/s')
ylabel('$y $',...
    'interpreter','latex',...
    'FontWeight','normal',...
    'FontSize',14,...
    'FontName','Times')
% daspect([1 1 1]);
txt = {'($a$)'};
text(1,0.5,txt,...
    'interpreter','latex',...
    'FontWeight','normal',...
    'FontSize',14,...
    'FontName','Times')
axis([0 max(x_1000km) 0 max(y_1000km)]);
% colorbar
set(gca,...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',14,...
    'FontName','Times')

%%
hh(2)=nexttile;
qp1=quiver(x_1000km(3:interval:end), y_1000km(3:interval:end), ...
    u_1(3:interval:end, 3:interval:end)',...
    v_1(3:interval:end, 3:interval:end)');
set(qp1,'AutoScale','on', 'AutoScaleFactor', 2,'color',[0 0 0])
legend('Scale: 2m/s')
xlabel('$x $'  ,'interpreter','latex',...
    'FontWeight','normal',...
    'FontSize',14,...
    'FontName','Times')
ylabel('$y $',...
    'interpreter','latex',...
    'FontWeight','normal',...
    'FontSize',14,...
    'FontName','Times')
txt = {'($b$)'};
text(1,0.5,txt,...
    'interpreter','latex',...
    'FontWeight','normal',...
    'FontSize',14,...
    'FontName','Times')
% daspect([1 1 1]);
axis([0 max(x_1000km) 0 max(y_1000km)]);
% colorbar
set(gca,...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',14,...
    'FontName','Times')
%%
sgtitle([height_title],...
    'interpreter','latex',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Times')
%%
tt.TileSpacing = 'compact';
tt.Padding = 'compact';