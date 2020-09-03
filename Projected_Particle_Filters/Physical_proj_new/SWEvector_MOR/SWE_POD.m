clear all;close all;clc
% load('SWE.mat')
load('SWE_run_4days.mat');
modeloutput= x_save;
[U,S,V] = svd(modeloutput,'econ');
N =size(modeloutput,1);
sig=diag(S);
figure(1)
plot(sig,'ko','Linewidth',(1.5)),grid on
xlabel('k')
ylabel('Singular value, \sigma_k')
% title('Standard plot of singular values')
%%
figure(2)
semilogy(sig,'-ok','LineWidth',1.5), grid on
xlabel('k')
ylabel('Semilogy of diag(S)')
% % cutoff =9e-9;
% r =max(find(sig>cutoff));% keep modes with sing.
% % figure
% % semilogy(diag(S),'-ok','LineWidth',1.5)
% % hold on,grid on
% % semilogy(diag(S(1:r,1:r)),'or','LineWidth',1.5)
% % plot([-20 N+20],[cutoff cutoff],'r--','LineWidth',2)
% % axis([0 6000 0 85550000])
% % xlabel('k')
% % ylabel('Semilogy of diag(S)')
% % 
% % 
% % figure
% % semilogy(diag(S),'-ok','LineWidth',1.5)
% % hold on,grid on
% % semilogy(diag(S(1:r,1:r)),'or','LineWidth',1.5)
% % plot([-20 N+20],[cutoff cutoff],'r--','LineWidth',2)
% % axis([0 800 0 100000000])
% % xlabel('k')
% % ylabel('Semilogy of diag(S)')
% 
% cdS =cumsum((diag(S)).^2)./sum((diag(S)).^2);% cumulative energy
% r90 =max(find(cdS>0.9));% find r that captures 90% energy
% X90 = U(:,1:r90)*S(1:r90,1:r90)*V(:,1:r90)';
% 
% figure
% semilogy(1-cdS,'-ok','LineWidth',1.5) 
% hold on,grid on
% semilogy(1-cdS(1:r90),'ob','LineWidth',1.5)
% % semilogy(1-cdS(1:r),'or','LineWidth',1.5)
% % set(gca,'XTick',[0 r90 300 600],'YTick',[0 0.8 0.99 1.0])
% % xlim([-10 600])
% % plot([r90 r90 0],[0 0.99 0.99],'b--','LineWidth',1.5)
% 
% 
% % r=20;
% % sig=diag(S);
% % figure(1)
% % plot(sig,'ko','Linewidth',(1.5)),grid on
% % xlabel('k')
% % ylabel('Singular value, \sigma_k')
% % % title('Standard plot of singular values')
% %% %%
% 
% % figure(2)
% % semilogy(sig,'bo','LineWidth',1.5), grid on
% % xlabel('k')
% % ylabel('Semilogy of diag(S)')
% % % hold off
% % % title('log plot of singular values')
% % 
% % % figure(3)
% % % semilogy(cdS,'ko','LineWidth',1),grid on
% % % xlabel('k')
% % % ylabel('Cumulative')
% % % hold on
% % nx =size(modeloutput,1); ny =size(modeloutput,2);
% % plotind = 2;
% % for r=[5 20 100]% truncation value
% % Xapprox = U(:,1:r)*S(1:r,1:r)*V(:,1:r)';% approximate image
% % subplot(2,2,plotind), 
% % plotind = plotind + 1;
% % imagesc(Xapprox),axis off
% % title(['r=',num2str(r,'%d'),', ',...
% % num2str(100*r*(nx+ny)/(nx*ny),'%2.2f'),'% storage']);
% % end
% % set(gcf,'Position',[100 100 550 400])
% % figure,subplot(1,2,1)
% % semilogy(diag(S),'k','LineWidth',1.2),grid on
% % xlabel('r')
% % ylabel('Singular value, \sigma_r')
% % xlim([-50 1550])
% % subplot(1,2,2)
% % plot(cumsum(diag(S))/sum(diag(S)),'k','LineWidth',1.2),grid on
% % xlabel('r')
% % ylabel('Cumulative Energy')
% % xlim([-50 1550])
% % ylim([0 1.1])
% % set(gcf,'Position',[100 100 550 240])
% % % Ur = U(:,1:r);
% % % Sr = S(1:r,1:r);
% % % Vr = V(:,1:r)'; % Truncate U,S,V using the rank r
% % % Xr = Ur*Sr*Vr; % Truncated matrix
% % % modeloutput_truncated=Xr;
% % % save('SWE_POD_r20.mat','modeloutput_truncated','x','y','H','dt','t_save')
% % 
% % % save('SWE_POD_r50_new.mat','modeloutput_truncated','x','y','H','dt','t_save')
