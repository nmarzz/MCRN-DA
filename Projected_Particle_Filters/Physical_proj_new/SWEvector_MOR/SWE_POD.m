clear all;close all;clc
% load('SWE.mat')
load('SWE_run_4days.mat');
modeloutput= x_save;
[U,S,V] = svd(modeloutput,'econ');
N =size(modeloutput,1);
sig=diag(S);
%%
% figure(1)
% plot(sig,'ko','Linewidth',(1.5)),grid on
% xlabel('k','fontsize',14,'interpreter','latex','FontName', 'Times New Roman','fontweight','bold')
% ylabel('Singular Values', 'fontsize',14,'interpreter','latex','FontName', 'Times New Roman','fontweight','bold')
%  set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
% % title('Standard plot of singular values')
TOLC=ptc12(12,'check');
figure(2)
semilogy(sig(1:350),'-ok','LineWidth',1.5), grid on, hold on
xlabel('i','fontsize',14,'interpreter','latex','FontName', 'Times New Roman','fontweight','bold')
ylabel('Singular Values $\sigma_{i}^{2}$','fontsize',14,'interpreter','latex','FontName', 'Times New Roman','fontweight','bold')
txt = ' $\longrightarrow R=10$' ;
% txt_1 = ' $\longrightarrow  r = 207$';
txt_2 =(  ' $\longrightarrow  R =30$');
txt_3 =(  ' $\longrightarrow  R =130$');
% txt_3 =(  ' $\nearrow $');
% txt_4 =(  ' $ R= 130$');
plot(10, sig(10),'o','Color', TOLC(10,:),'LineWidth',1.5)
% plot(207, sig(207),'o','Color', TOLC(4,:),'LineWidth',1.5)
plot(30, sig(30),'o','Color', TOLC(10,:),'LineWidth',1.5)
plot(130, sig(130),'o','Color', TOLC(10,:),'LineWidth',1.5)
h=text(10,sig(10),txt);
% h_1=text(207,sig(207),txt_1);
h_2=text(30,sig(30),txt_2);
 h_3=text(130,sig(130),txt_3);
%  h_4=text(130,sig(130),txt_4);
set([h,h_2,h_3], 'Color', TOLC(10,:),'interpreter','latex','fontsize',11,'FontName', 'Times New Roman','fontweight','bold')
% ylim([1.000000000000000e-16, 10^09]);
xlim([0,350]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
hold off
% 
% figure(3)
% semilogy(sig,'-ok','LineWidth',1.5), grid on
% xlabel('k','fontsize',14,'interpreter','latex','FontName', 'Times New Roman','fontweight','bold')
% ylabel('Singular Values','fontsize',14,'interpreter','latex','FontName', 'Times New Roman','fontweight','bold')
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
% %  title('log plot of singular values')
%% 
cdS =cumsum((sig.^2)./sum(sig.^2));% cumulative 
error=1-cdS;
figure(4)
semilogy(error,'k-o','LineWidth',1),grid on,hold on;
% txt = (' Error=$$\sqrt{\sum_{i=R+1}^{K} \sigma_{i}^{2}}$$') ;
% h_1=text(150,0.0000001,txt);
txt_1= ' $\longrightarrow R=10$' ;
txt_2 =(  ' $\longrightarrow  R =30$');
txt_3 =(  ' $\longrightarrow  R =130$');
plot(10, error(10),'o','Color', TOLC(10,:),'LineWidth',1.5)
plot(30, error(30),'o','Color', TOLC(10,:),'LineWidth',1.5)
plot(130, error(130),'o','Color', TOLC(10,:),'LineWidth',1.5)
h_1=text(10,error(10),txt_1);
h_2=text(30,error(30),txt_2);
h_3=text(130,error(130),txt_3);
set([h_1,h_2,h_3], 'Color', TOLC(10,:),'interpreter','latex','fontsize',11,'FontName', 'Times New Roman','fontweight','bold')
xlabel('i','fontsize',14,'interpreter','latex','FontName', 'Times New Roman','fontweight','bold')
ylabel('Error','fontsize',14,'interpreter','latex','FontName', 'Times New Roman','fontweight','bold')
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)

xlims = get(gca,'XLim');
plot(xlims, [error(10) error(10)+10], '-k')
plot(xlims, [10 10], '-r')




% r=30;
% Ur = U(:,1:r);
% Sr = S(1:r,1:r);
% Vr = V(:,1:r)'; % Truncate U,S,V using the rank r
% Xr = Ur*Sr*Vr; % Truncated matrix
% modeloutput_truncated=Xr;
% % % % save('SWE_POD_r20.mat','modeloutput_truncated','x','y','H','dt','t_save')
% save('SWE_POD_r30_new.mat','modeloutput_truncated','x','y','H','dt','t_save')
