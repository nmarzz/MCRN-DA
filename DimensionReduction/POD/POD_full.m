clc;clear all;close all;
load('rotor_oscillator_vf.mat')
whos
%The data have:
%Ux-Horizontal velocity components
%Uy-Vertical velocity components
%x-dim,y-dim for each dt-time step and t is the time

%% Step1,build a matrix such that time x-dir and space y-dir 
L=length(t);%number of snapshots 
m=length(y);%x-dimension
n=length(x);%y-dimension 
k=m*n;%each snapshot has x*y-dimensions
Ux_reshaped = reshape( Ux, [], L);
Uy_reshaped = reshape( Uy, [], L);
X=[Ux_reshaped;Uy_reshaped];
%% POD
[U,S,V]=svd(X);
Tol=0.99;
sig=diag(S);
cdS =cumsum(sig.^2)./sum(sig.^2);% cumulative energy
r =find(cdS>Tol, 1 ); 

figure(1)
plot(sig,'ko','Linewidth',(1.5)),grid on
xlabel('i')
ylabel('Singular value, \sigma_i')
title('Standard plot of singular values')

figure(2)
semilogy(diag(S),'bo','LineWidth',1.5), grid on
xlabel('i')
ylabel('Semilogy of diag(S)')
hold off
title('log plot of singular values')

figure(3)
plot(cdS,'ko','LineWidth',1.2),grid on
xlabel('i')
ylabel('Cumulative Energy')
%% Matrix Truncation
U1=U(:,1:r);S1=S(1:r,1:r);V1=V(:,1:r)';%regenerate U,S,V using the rank r
Xr=U1*S1*V1; %Truncated matrix
Uxr_reshaped=Xr(1:k,:);
Uyr_reshaped=Xr(k+1:end,:); 

Uxr_H_Velocity_truncated=reshape(Uxr_reshaped,m,n,L);
Uyr_V_Velocity_truncated=reshape(Uyr_reshaped,m,n,L);
%% visulazation
slice=5;
figure(4)
waterfall(x,y,Ux(:,:,slice).^2 + Uy(:,:,slice).^2),set(gca,'xlim',[-3.5 4],'Zlim',[0 8])
title ('Representation of the speed surface of X' )
figure(5)
waterfall(x,y,Uxr_H_Velocity_truncated(:,:,slice).^2 + Uyr_V_Velocity_truncated(:,:,slice).^2), set(gca,'xlim',[-3.5 4],'Zlim',[0 8])
title ('Representation of the speed surface of X_r' )

figure(6)
t = tiledlayout(2,2,'TileSpacing','Compact');
nexttile
h1=quiver( x, y, Ux(:,:,slice),Uy(:,:,slice)); 
set(h1,'AutoScale','on', 'AutoScaleFactor', 1.5,'color',[0 0 1])
set(gca,'xlim',[-0.8 0.8],'Ylim',[0 1])
title('X')
%
nexttile
h1=quiver( x, y, Uxr_H_Velocity_truncated(:,:,slice),Uyr_V_Velocity_truncated(:,:,slice),'color',[1 0 1]); 
set(h1,'AutoScale','on', 'AutoScaleFactor', 1.5)
set(gca,'xlim',[-0.8 0.8],'Ylim',[0 1])
title('X_r')
% 
nexttile([1 2])
h1=quiver( x, y, Ux(:,:,slice),Uy(:,:,slice),'color',[0 0 1]); 
hold on
h2=quiver( x, y, Uxr_H_Velocity_truncated(:,:,slice),Uyr_V_Velocity_truncated(:,:,slice),'color',[1 0 1]); 
set([h1 h2],'AutoScale','on', 'AutoScaleFactor', 1.5)
set(gca,'xlim',[-0.8 0.8],'Ylim',[0 1])
title('Comparison')
title(t,'\fontsize{16}Velocity vector field of X & X_r ')
%% save the truncated data as mat file
save('velocity','Uxr_H_Velocity_truncated','Uyr_V_Velocity_truncated','x','y','t','dt')
