clc;clear all;close all;
load('rotor_oscillator_vf.mat')
whos
%The data have:
%Ux-Horizontal velocity components
%Uy-Vertical velocity components
%x-dim,y-dim for each dt-time step and t is the time

%% Step1, build a matrix such that time x-dir and space y-dir
L=length(t);%number of snapshots
m=length(y);%x-dimension
n=length(x);%y-dimension
N_gridpoints=m*n;%each snapshot has x*y-dimensions
Ux_reshaped = reshape( Ux, [], L);
Uy_reshaped = reshape( Uy, [], L);
X=[Ux_reshaped;Uy_reshaped];
%% POD
[U,S,V]=svd(X,'econ');
Tol=0.99;
sig=diag(S);
cdS =cumsum(sig.^2)./sum(sig.^2);% cumulative
r =find(cdS>Tol, 1 );

figure(1)
plot(sig,'ko','Linewidth',(1.5)),grid on
xlabel('i')
ylabel('Singular value, \sigma_i')
title('Standard plot of singular values')
%
figure(2)
semilogy(sig,'bo','LineWidth',1.5), grid on
xlabel('i')
ylabel('Semilogy of diag(S)')
title('log plot of singular values')
% 
figure(3)
plot(cdS,'ko','LineWidth',1.2),grid on
xlabel('i')
ylabel('Cumulative')


%% Matrix Truncation
U_r=U(:,1:r);S_r=S(1:r,1:r);V_r=transpose(V(:,1:r));%regenerate U,S,V using the rank r
Xr=U_r*S_r*V_r; %Truncated matrix

Uxr_reshaped=Xr(1:N_gridpoints,:);
Uyr_reshaped=Xr(N_gridpoints+1:end,:);

Uxr_H=reshape(Uxr_reshaped,m,n,L);%Horizontal velocity truncated
Uyr_V=reshape(Uyr_reshaped,m,n,L);%Vertical velocity truncated
%% visulazation
slice=4;
figure(5)
%plot original velocity
waterfall(x,y,Ux(:,:,slice).^2 + Uy(:,:,slice).^2),set(gca,'xlim',[-3.5 4],'Zlim',[0 8])
title ('Representation of the speed surface of X' )
figure(6)
%plot truncated velocity
waterfall(x,y,Uxr_H(:,:,slice).^2 + Uyr_V(:,:,slice).^2), set(gca,'xlim',[-3.5 4],'Zlim',[0 8])
title ('Representation of the speed surface of X_r' ) 
figure(7)
surfl(x,y,Ux(:,:,slice));shading interp;
title ('spatial (velocity) field Ux' )
%% Comparison
figure(9)
t = tiledlayout(2,1,'TileSpacing','Compact');
nexttile
%plot original velocity
h1=quiver( x, y, Ux(:,:,slice),Uy(:,:,slice)); 
set(h1,'AutoScale','on', 'AutoScaleFactor', 1.5,'color',[0 0 1])
set(gca,'xlim',[-0.8 0.8],'Ylim',[0 1])
title('X')
nexttile
%plot truncated velocity
h2=quiver( x, y, Uxr_H(:,:,slice),Uyr_V(:,:,slice),'color',[1 0 1]);
set(h2,'AutoScale','on', 'AutoScaleFactor', 1.5)
set(gca,'xlim',[-0.8 0.8],'Ylim',[0 1])
title('X_r')
title(t,'\fontsize{16}Velocity vector field of X & X_r ')
%% Error between X and X_r
Error=X-Xr;%difference matrix between X and X_r
Ux_E_reshaped=Error(1:N_gridpoints,:);
Uy_E_reshaped=Error(N_gridpoints+1:end,:);
Ux_E=reshape(Ux_E_reshaped,m,n,L);%Horizontal velocity of the difference matrix
Uy_E=reshape(Uy_E_reshaped,m,n,L);%Vertical velocity of the difference matrix
figure(10)
f=tiledlayout(2,1,'TileSpacing','Compact');
nexttile
%plot original velocity
h1=quiver( x, y, Ux(:,:,slice),Uy(:,:,slice));
set(h1,'AutoScale','on', 'AutoScaleFactor', 1.5,'color',[0 0 1])
hold on
%plot truncated velocity
h2=quiver( x, y, Uxr_H(:,:,slice),Uyr_V(:,:,slice));
set(gca,'xlim',[-0.8 0.8],'Ylim',[0 1])
set(h2,'AutoScale','on', 'AutoScaleFactor', 1.5,'color',[1 0 1])
title('Comparison')
nexttile
%plot the error difference
hh=quiver( x, y,Ux_E(:,:,slice),Uy_E(:,:,slice));
set(hh,'AutoScale','on', 'AutoScaleFactor', 1.5,'color',[0 0 1])
set(gca,'xlim',[-0.8 0.8],'Ylim',[0 1])
title('Difference')
title(f,'\fontsize{16}Comparison between velocity vector field of X & X_r ')
%% different spatial mode of original matrix X
%U is a spatial structure from svd
Ux_m=U(1:N_gridpoints,:);
Uy_m=U(N_gridpoints+1:end,:);
%% Spatial modes of the horizontal velocity
figure(11)
ff=tiledlayout(2,2,'TileSpacing','Compact');
%mode1
nexttile
plot(Ux_m(:,1),'k','Linewidth',(1.5)); 
set(gca,'Fontsize',(8))
title('mode 1')
%mode2
nexttile
plot(Ux_m(:,2),'--m','Linewidth',(1.5));
set(gca,'Fontsize',(8))
title('mode 2')
%mode1,2,3
nexttile([1 2])
plot(Ux_m(:,1),'k','Linewidth',(1.5)); 
hold on
plot(Ux_m(:,2),'--m','Linewidth',(1.5));%mode2
plot(Ux_m(:,3),'.c','Linewidth',(1.5));%mode3
legend('mode 1','mode 2','mode 3','Location','SouthEast')
set(gca,'Fontsize',(8))
hold off
title (ff,'Spatial modes of the horizontal velocity' )
%% Spatial modes of the Vertical velocity
figure(12)
fff=tiledlayout(2,2,'TileSpacing','Compact');
%mode1
nexttile
plot(Uy_m(:,1),'k','Linewidth',(1.5)); 
set(gca,'Fontsize',(8))
title('mode 1')
%mode2
nexttile
plot(Uy_m(:,2),'--m','Linewidth',(1.5));
set(gca,'Fontsize',(8))
title('mode 2')
%mode1,2,3
nexttile([1 2])
plot(Uy_m(:,1),'k','Linewidth',(1.5)); 
hold on
plot(Uy_m(:,2),'--m','Linewidth',(1.5));%mode2
plot(Uy_m(:,3),'.c','Linewidth',(1.5));%mode3
legend('mode 1','mode 2','mode 3','Location','SouthEast')
set(gca,'Fontsize',(8))
hold off
title (fff,'Spatial modes of the vertical velocity' )

%% time coefficient mode of matrix X
figure(13)
plot(t,V(:,1),'k',t,V(:,2),'-c','Linewidth',(2))
set(gca,'Fontsize',(8))
legend('mode 1','mode 2','Location','SouthEast')

