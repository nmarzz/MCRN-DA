clear all;close all;clc
% load('SWE.mat')
load('SWE_run_new.mat');
modeloutput= x_save;
[U,S,V] = svd(modeloutput,'econ');
r=20;
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
% save('SWE_POD_r20.mat','modeloutput_truncated','x','y','H','dt','t_save')

save('SWE_POD_r20_new.mat','modeloutput_truncated','x','y','H','dt','t_save')
