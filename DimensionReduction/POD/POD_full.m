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
Ux1 = reshape( Ux, [], L);
Uy1 = reshape( Uy, [], L);
X=[Ux1;Uy1];
%% POD
[U,S,V]=svd(X);
Tol=0.;
sig=diag(S);
cdS =cumsum(sig.^2)./sum(sig.^2);% cumulative energy
r =find(cdS>Tol, 1 ); 

% figure(1)
% plot(sig,'ko','Linewidth',(1.5)),grid on
% xlabel('k')
% ylabel('Singular value, \sigma_k')
% title('Standard plot of singular values')
% 
% figure(2)
% semilogy(diag(S),'bo','LineWidth',1.5), grid on
% xlabel('k')
% ylabel('Semilogy of diag(S)')
% hold off
% title('log plot of singular values')
% 
% figure(3)
% plot(cdS,'ko','LineWidth',1.2),grid on
% xlabel('k')
% ylabel('Cumulative Energy')
%% Matrix Truncation
U1=U(:,1:r);S1=S(1:r,1:r);V1=V(:,1:r)';%regenerate U,S,V using the rank r
X1=U1*S1*V1; %Truncated matrix
Ux2=X1(1:k,:);
Uy2=X1(k+1:end,:); 

Ux3=reshape(Ux2,m,n,L);
Uy3=reshape(Uy2,m,n,L);
%% visulazation

figure(4)
surfl(x,y,Ux(:,:,1).^2+Uy(:,:,1).^2);shading interp;
title ('spatial (velocity) field Ux' )
figure(5)
surfl(x,y,Ux3(:,:,end));shading interp;
title ('spatial (velocity) field Ux_r' )
%
figure(6)
surfl(x,y,Ux(:,:,end));shading interp;
title ('Time-evolution coefficient Uy' )
figure(7)
surfl(x,y,Ux3(:,:,end));shading interp;
title ('Time-evolution coefficient Uy_r' )
figure(8)
waterfall(x,y,Ux(:,:,1).^2 + Uy(:,:,1).^2),set(gca,'xlim',[-3.5 4],'Zlim',[0 8])
title ('Representation of the speed surface of X' )
figure(9)
waterfall(x,y,Ux3(:,:,1).^2 + Uy3(:,:,1).^2), set(gca,'xlim',[-3.5 4],'Zlim',[0 8])
title ('Representation of the speed surface of X_r' )


