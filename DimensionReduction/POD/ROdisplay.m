clear all;close all;clc;
%Preprocessing
% load('velocity.mat');
% Ux=Ux3;
% Uy=Uy3;

% load the velocity field
load rotor_oscillator_vf.mat

% list the variables
whos

% restrict the domain to [-2,2] x [-1,1]
selx = abs(x) <= 2;
sely = abs(y) <= 1;
x = x(selx);
y = y(sely);
Ux = Ux(:,selx,:);
Ux = Ux(sely,:,:);
Uy = Uy(:,selx,:);
Uy = Uy(sely,:,:);

%% feel free to process this restricted field Ux Uy using modal decomposition
%% the field outside of the window is not interesting to us

%%
% Calculate some interesting physical quantities

% compute pointwise speed and direction vector field
% (all vectors are unit length, but orientation is the same as for velocity)
S = hypot(Ux,Uy);
Dirx = Ux./S;
Diry = Uy./S;

% compute vorticity (up to a constant )
[Uxx,Uxy] = gradient(Ux);
[Uyx,Uyy] = gradient(Uy);
Vort = Uyx-Uxy;
Div = Uxx + Uyy;
logDiv = log10(abs(Div)); % we expect small divergence

% divergent colormap (for visualizing vorticity)
% see: https://blogs.mathworks.com/steve/2015/01/20/divergent-colormaps/
div_colors = [(winter(64)); flipud(autumn(64))];

disp('Visualizing the data');

% figure(1);
% filename ='frames.gif';
for ti = 1:1073
    
    if ti == 1 % in first step draw the arrows
        figure(1);
        [~,hs] = contourf(x,y, log10(S(:,:,ti)), 50);
        
        hold on;
        hq = quiver(x,y, ...
            Dirx(:,:,ti), ...
            Diry(:,:,ti) ); hold off;
        
        % red arrows that are cut off by screen
        hq.Color = 'r'; hq.Clipping = 'on';
        
        fig1=gca;
        xlabel(fig1,'X');ylabel(fig1,'Y');
        
        % fixed axes
        axis equal;axis tight;
        axis(fig1,[-2,2,-1,1]);
        xlim(fig1,'manual');
        ylim(fig1,'manual');
        
        set(fig1,'color',repmat(0.7,[1,3])); %gray background
        colormap(parula);
        hsc = colorbar; %title(hsc,'log_{10} Speed');
        
    else % in following steps, just reposition arrows
        hs.ZData = S(:,:,ti);
        hq.UData = Dirx(:,:,ti);
        hq.VData = Diry(:,:,ti);
    end
    title(fig1,sprintf('Velocity snapshot= %.f',ti));
%     title(fig1,sprintf('Velocity/speed at t = %.2f',t(ti)));
    
%     if ti == 1 % in first step create color plot
%         figure(2);
%         [~,hv] = contourf(x,y, Vort(:,:,ti)); shading flat;
%         
%         fig2=gca;
%         set(fig2,'color',repmat(0.7,[1,3])); %gray background
%         xlabel('X');ylabel('Y');
%         axis equal; axis tight; colorbar;
%         caxis([-2,2]);
%         colormap(div_colors);
%     else % in following steps, just reposition arrows
%         hv.ZData = Vort(:,:,ti);
%     end
%     title(fig2,sprintf('Vorticity at t = %.2f',t(ti)));
%     
%     if ti == 1 % in first step create color plot
%         figure(3);
%         [~,hdiv] = contourf(x,y, logDiv(:,:,ti)); shading flat;
%         
%         fig3=gca;
%         set(fig3,'color',repmat(0.7,[1,3])); %gray background
%         xlabel('X');ylabel('Y');
%         axis equal; axis tight; colorbar;
%         hdc = colorbar; title(hdc,'log_{10} Divergence');
%     else % in following steps, just reposition arrows
%         
%         hdiv.ZData = Div(:,:,ti);
%     end
%     title(fig3,sprintf('Divergence at t = %.2f',t(ti)));
%     pause(0.0001)
%     frame = getframe(gcf);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     % Write to the GIF File
%     if ti == 1
%         imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',1);
%     else
%         imwrite(imind,cm,filename,'gif','WriteMode','append');
%     end
    
    
end
