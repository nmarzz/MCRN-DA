% This script animates the height field and the vorticity produced by
% a shallow water model. It should be called only after shallow_water_model
% has been run.
clc;clear all;close all;
load('SWE.mat')
modeloutput=x_save; % To plot the original SWE
% load('SWE_truncated.mat')
% modeloutput=modeloutput_truncated;
L=length(t_save);%number of snapshots
m=length(y);%y-dimension
n=length(x);%x-dimension
N_gridpoints=m*n;%each snapshot has x*y-dimensions
u=modeloutput(1:N_gridpoints,:);v=modeloutput(N_gridpoints+1:N_gridpoints*2,:);h=modeloutput((N_gridpoints*2)+1:end,:);
u_save=reshape(u,n,m,L);v_save=reshape(v,n,m,L);h_save=reshape(h,n,m,L);

%%
% Axis units are thousands of kilometers (x and y are in metres)
x_1000km = x.*1e-6;
y_1000km = y.*1e-6;

% Set colormap to have 64 entries
ncol=64;
colormap(jet(ncol));

% Interval between arrows in the velocity vector plot
interval = 6;

% Set this to "true" to save each frame as a png file
plot_frames = false;

% Decide whether to show height in metres or km
if mean(plot_height_range) > 1000
    height_scale = 0.001;
    height_title = 'Height (km)';
else
    height_scale = 1;
    height_title = 'Height (m)';
end

disp(['Maximum orography height = ' num2str(max(H(:))) ' m']);

% Loop through the frames of the animation
% it = 1:noutput
snapshot=15;
for it =snapshot
    clf
    
    % Extract the height and velocity components for this frame
    h = squeeze(h_save(:,:,it));
    u = squeeze(u_save(:,:,it));
    v = squeeze(v_save(:,:,it));
    
    % Plot the height field
    handle = image(x_1000km, y_1000km, (h'+H').*height_scale);
    set(handle,'CDataMapping','scaled');
    set(gca,'ydir','normal');
    caxis(plot_height_range.*height_scale);
    
    % Plot the orography as black contours every 1000 m
    hold on
    warning off
    contour(x_1000km, y_1000km, H',[1:1000:8001],'k');
    warning on
    
    % Plot the velocity vectors
    quiver(x_1000km(3:interval:end), y_1000km(3:interval:end), ...
        u(3:interval:end, 3:interval:end)',...
        v(3:interval:end, 3:interval:end)','k');
    
    % Write the axis labels, title and time of the frame
    xlabel('X distance (1000s of km)');
    ylabel('Y distance (1000s of km)');
    title(['\bf' height_title]);
    text(0, max(y_1000km), ['Time=' num2str(t_save(it)./3600) 'hours'],...
        'verticalalignment','bottom','fontsize',12);
    % Set other axes properties and plot a colorbar
%     daspect([2 .5 2]);
    axis([0 max(x_1000km) 0 max(y_1000km)]);
    colorbar
end
