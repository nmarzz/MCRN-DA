clear all;close all;clc
load('SWE.mat')
% whos
modeloutput=x_save(:,1440:end);%snapshot matrix
t_save=t_save(:,1440:end);
L=length(t_save);%number of snapshots
m=length(y);%x-dimension
n=length(x);%y-dimension
N_gridpoints=m*n;%each snapshot has x*y-dimensions
[U,S,V] = svd(modeloutput,'econ');
r=2;
sig=diag(S);
% figure(1)
% plot(sig,'ko','Linewidth',(1.5)),grid on
% xlabel('k')
% ylabel('Singular value, \sigma_k')
% title('Standard plot of singular values')

% figure(2)
% semilogy(diag(S),'bo','LineWidth',1.5), grid on
% xlabel('k') 
% ylabel('Semilogy of diag(S)')
% hold off
% title('log plot of singular values')

cdS =cumsum(sig.^2)./sum(sig.^2);% cumulative
% figure(3)
% plot(cdS,'ko','LineWidth',1),grid on
% xlabel('k')
% ylabel('Cumulative')

Ur = U(:,1:r);
Sr = S(1:r,1:r);
Vr = V(:,1:r)'; % Truncate U,S,V using the rank r
Xr = Ur*Sr*Vr; % Truncated matrix
modeloutput_truncated=Xr;
% save('SWE_truncated','modeloutput_truncated','x','y','H','dt','t_save','plot_height_range')

