% This script animates the height field and the vorticity produced by
% a shallow water model. It should be called only after shallow_water_model
% has been run.
clear all; close all; clc;
% load('SWE.mat')
% % clear all; close all; clc;
% load('SWE_run_4days.mat');
% modeloutput= x_save;
%% POD
% load('SWE_POD_r2.mat')
% load('SWE_POD_r10.mat')
load('SWE_POD_r30_new.mat')
modeloutput=modeloutput_truncated;
%% DMD
% load('SWE_DMD_r2.mat')
% load('SWE_DMD_r10.mat')
% load('SWE_DMD_r20.mat')
% modeloutput=X_dmd;

L=length(t_save);%number of snapshots
m=length(y);%y-dimension
n=length(x);%x-dimension
N_gridpoints=m*n;%each snapshot has x*y-dimensions

u=modeloutput(1:N_gridpoints,:);v=modeloutput(N_gridpoints+1:N_gridpoints*2,:);h=modeloutput((N_gridpoints*2)+1:end,:);
u_save=reshape(u,n,m,L);v_save=reshape(v,n,m,L);h_save=reshape(h,n,m,L);

%%
% This script animates the height field and the vorticity produced by
% a shallow water model. It should be called only after shallow_water_model
% has been run.
% Set the size of the figure
set(gcf,'units','inches');
pos=get(gcf,'position');
pos([3 4]) = [10.5 5];
set(gcf,'position',pos)
plot_height_range = [9500 10500];
% Set other figure properties and draw now
set(gcf,'defaultaxesfontsize',12,...
    'paperpositionmode','auto','color','w');
drawnow

% Axis units are thousands of kilometers (x and y are in metres)
x_1000km = x.*1e-6;
y_1000km = y.*1e-6;

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
    height_title = 'POD(r=30)';
else
    height_scale = 1;
    height_title = 'POD(r=30)';
end

disp(['Maximum orography height = ' num2str(max(H(:))) ' m']);

% Loop through the frames of the animation
filename = 'testnew51.gif';
for it =1:60: L
    clf
    
    % Extract the height and velocity components for this frame
    h = squeeze(h_save(:,:,it));
    u = squeeze(u_save(:,:,it));
    v = squeeze(v_save(:,:,it));
    
    % First plot the height field
    %   subplot(2,1,1);
    
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
%     xlabel('X distance (1000s of km)');
%     ylabel('Y distance (1000s of km)');
    title(['\bf' height_title]);
%     text(0, max(y_1000km), ['Time = ' num2str(t_save(it)./3600) ' hours'],...
%         'verticalalignment','bottom','fontsize',12);
    
    % Set other axes properties and plot a colorbar
    daspect([1 0.5 1]);
    axis([0 max(x_1000km) 0 max(y_1000km)]);
    colorbar
    
    %   % Compute the vorticity
    %   vorticity = zeros(size(u));
    %   vorticity(2:end-1,2:end-1) = (1/dy).*(u(2:end-1,1:end-2)-u(2:end-1,3:end)) ...
    %      + (1/dx).*(v(3:end,2:end-1)-v(1:end-2,2:end-1));
    
    % Now plot the vorticity
    %   subplot(2,1,2);
    %   handle = image(x_1000km, y_1000km, vorticity');
    %   set(handle,'CDataMapping','scaled');
    %   set(gca,'ydir','normal');
    %   caxis([-3 3].*1e-4);
    
    %     % Axis labels and title
    %     xlabel('X distance (1000s of km)');
    %     ylabel('Y distance (1000s of km)');
    %     title('\bfRelative vorticity (s^{-1})');
    
    %     % Other axes properties and plot a colorbar
    %     daspect([1 1 1]);
    %     axis([0 max(x_1000km) 0 max(y_1000km)]);
    %     colorbar
    
    % Now render this frame
    warning off
    drawnow
    warning on
    
    % To make an animation we can save the frames as a
    % sequence of images
%         if plot_frames
%             imwrite(frame2im(getframe(gcf)),...
%                 ['frame' num2str(it,'%03d') '.png']);
%         end
    pause(1)
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File
    if it == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',1);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
end

% % Axis units are thousands of kilometers (x and y are in metres)
% % Set the size of the figure
% set(gcf,'units','inches');
% pos=get(gcf,'position');
% pos([3 4]) = [10.5 5];
% set(gcf,'position',pos)
% plot_height_range = [9500 10500];
% % Set other figure properties and draw now
% % set(gcf,'defaultaxesfontsize',10,...
% %     'paperpositionmode','auto','color','w');
% drawnow
%
% % Axis units are thousands of kilometers (x and y are in metres)
% x_1000km = x.*1e-6;
% y_1000km = y.*1e-6;
%
% % Set colormap to have 64 entries
% ncol=64;
% colormap(jet(ncol));
%
% % Interval between arrows in the velocity vector plot
% interval = 10;
%
% % Set this to "true" to save each frame as a png file
% plot_frames = false;
%
% % Decide whether to show height in metres or km
% if mean(plot_height_range) > 1000
%     height_scale = 0.001;
%     height_title = 'Height (km)';
% else
%     height_scale = 1;
%     height_title = 'Height (m)';
% end
%
% disp(['Maximum orography height = ' num2str(max(H(:))) ' m']);
% %% Loop through the frames of the animation
% % prevFontName = get(0,'defaultAxesFontName');
% % set(0,'defaultAxesFontName','Arial');
% % resolutionScaling = 96/get(0,'ScreenPixelsPerInch');
% % fontSize = 16*resolutionScaling;
% % lineWidth = 2*resolutionScaling;
% % markerSize = 10*resolutionScaling;
% % for it = 10
% %     clf
% %     figure
% %
% %     % Extract the height and velocity components for this frame
% %     h = squeeze(h_save(:,:,it));
% %     u = squeeze(u_save(:,:,it));
% %     v = squeeze(v_save(:,:,it));
% %     % Plot the height field
% %     handle = image(x_1000km, y_1000km, (real(h)'+H').*height_scale);
% %     set(handle,'CDataMapping','scaled');
% %     caxis(plot_height_range.*height_scale);
% %
% %     % Plot the orography as black contours every 1000 m
% %     hold on
% %     warning off
% %     contour(x_1000km, y_1000km, H',[1:1000:8001],'k','LineWidth',lineWidth);
% %     warning on
% %
% %     % Plot the velocity vectors
% %     quiver(x_1000km(3:interval:end), y_1000km(3:interval:end), ...
% %         real(u(3:interval:end, 3:interval:end))',...
% %         real(v(3:interval:end, 3:interval:end))','k','LineWidth',lineWidth);
% %     set(gca,'Fontsize',fontSize,'LineWidth',lineWidth)
% %
% %     % Write the axis labels, title and time of the frame
% %     xlabel('X distance (1000s of km)','Fontsize',fontSize);
% %     ylabel('Y distance (1000s of km)','Fontsize',fontSize);
% %     title([height_title]);
% %     text(0, max(y_1000km), ['Time = ' num2str(t_save(it)./3600) ' hours'],...
% %         'verticalalignment','bottom','Fontsize',fontSize);
% %     set(0,'defaultAxesFontName',prevFontName);
% %     % Set other axes properties and plot a colorbar
% %     daspect([1 1 1]);
% %     axis([0 max(x_1000km) 0 max(y_1000km)]);
% %     colorbar
% %     fig2svg('try2d.svg'); %  ,'',1);
% %
% %
% %     % Now render this frame
% %     warning off
% %     drawnow
% %     warning on
% %
% %     % To make an animation we can save the frames as a
% %     % sequence of images
% %     if plot_frames
% %         imwrite(frame2im(getframe(gcf)),...
% %             ['frame' num2str(it,'%03d') '.png']);
% %     end
% % end

