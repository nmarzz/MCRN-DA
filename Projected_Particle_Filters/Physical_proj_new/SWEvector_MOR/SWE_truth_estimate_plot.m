clear all; close all; clc;
load('SWE_run_4days.mat');
% modeloutput= x_save;
load('SWE_POD.mat')
modeloutput= ess;
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
    height_title = 'Height ($h$), Truth';
else
    height_scale = 1;
    height_title = 'Height ($h$), Truth';
end
it=60*24;
% Extract the height and velocity components for this frame
h = squeeze(h_save(:,:,it));
u = squeeze(u_save(:,:,it));
v = squeeze(v_save(:,:,it));
% Plot the height field
it_1=60*48;
% Extract the height and velocity components for this frame
h_1= squeeze(h_save(:,:,it_1));
u_1 = squeeze(u_save(:,:,it_1));
v_1 = squeeze(v_save(:,:,it_1));
it_2=60*72;
% Extract the height and velocity components for this frame
h_2= squeeze(h_save(:,:,it_2));
u_2 = squeeze(u_save(:,:,it_2));
v_2 = squeeze(v_save(:,:,it_2));
%   warning on
figure(1)
tt = tiledlayout(3,1);
hh(1)=nexttile;
handle = image(x_1000km, y_1000km, (h'+H').*height_scale);
set(handle,'CDataMapping','scaled');
set(gca,'ydir','normal');
caxis(plot_height_range.*height_scale);
% xlabel('x '  ,'FontUnits','points',...
%     'FontWeight','normal',...
%     'FontSize',14,...
%     'FontName','Times');
% ylabel('$y $',...
%     'interpreter','latex',...
%     'FontWeight','normal',...
%     'FontSize',14,...
%     'FontName','Times')
daspect([1 1 1]);
axis([0 max(x_1000km) 0 max(y_1000km)]);
% colorbar
set(gca,...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',14,...
    'FontName','Times')
%%
hh(2)=nexttile;
handle = image(x_1000km, y_1000km, (h_1'+H').*height_scale);
set(handle,'CDataMapping','scaled');
set(gca,'ydir','normal');
caxis(plot_height_range.*height_scale);
% xlabel('x '  ,'FontUnits','points',...
%     'FontWeight','normal',...
%     'FontSize',14,...
%     'FontName','Times');
ylabel('$y $',...
    'interpreter','latex',...
    'FontWeight','normal',...
    'FontSize',14,...
    'FontName','Times')
daspect([1 1 1]);
axis([0 max(x_1000km) 0 max(y_1000km)]);
% colorbar
set(gca,...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',14,...
    'FontName','Times')
%%
hh(3)=nexttile;
handle = image(x_1000km, y_1000km, (h_2'+H').*height_scale);
set(handle,'CDataMapping','scaled');
set(gca,'ydir','normal');
caxis(plot_height_range.*height_scale);
xlabel('$x $' ,...
    'interpreter','latex',...
    'FontWeight','normal',...
    'FontSize',14,...
    'FontName','Times')

% ylabel('$y $',...
%     'interpreter','latex',...
%     'FontWeight','normal',...
%     'FontSize',14,...
%     'FontName','Times')

daspect([1 1 1]);
axis([0 max(x_1000km) 0 max(y_1000km)]);
% colorbar
set(gca,...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',14,...
    'FontName','Times')

sgtitle([height_title],...
    'interpreter','latex',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Times')
%%
set(hh, 'Colormap', colormap(cmocean('curl')))
cbh = colorbar(hh(2)); 
% Reposition to figure's left edge, centered vertically
cbh.Position(1) = .95-cbh.Position(3);
cbh.Position(2) = 0.5-cbh.Position(4)/2;
% decrease horizontal extent of subplots to 92% of their current width
set(hh, {'Position'}, mat2cell(vertcat(hh.Position) .* [1 1 .92, 1], ones(size(hh(:))),4))
% tt.TileSpacing = 'compact';
% tt.Padding = 'compact';