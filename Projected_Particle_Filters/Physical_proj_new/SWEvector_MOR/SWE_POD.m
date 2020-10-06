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
txt = ' $\longrightarrow R=60$' ;
txt_1 = ' $\longrightarrow  R = 30$';
% txt_2 =(  ' $\longrightarrow  R =305$');
% txt_3 =(  ' $\longrightarrow  R =130$');
txt_3 =(  ' $\nearrow $');
txt_4 =(  ' $ R= 305$');
plot(60, sig(60),'o','Color', TOLC(10,:),'LineWidth',1.5)
plot(30, sig(30),'o','Color', TOLC(10,:),'LineWidth',1.5)
plot(305, sig(305),'o','Color', TOLC(10,:),'LineWidth',1.5)
% plot(130, sig(130),'o','Color', TOLC(10,:),'LineWidth',1.5)
h=text(60,sig(60),txt);
h_1=text(30,sig(30),txt_1);
% h_2=text(300,sig(300),txt_2);
 h_3=text(302,sig(280),txt_3);
 h_4=text(305,sig(230),txt_4);
set([h,h_1,h_3,h_4], 'Color', TOLC(10,:),'interpreter','latex','fontsize',11,'FontName', 'Times New Roman','fontweight','bold')
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
txt_1= ' $\longrightarrow R=30$' ;
txt_2 =(  ' $\longrightarrow  R =60$');
% txt_3 =(  ' $\longrightarrow  R =305$');
txt_3 =(  ' $\nearrow $');
txt_4 =(  ' $ R= 305$');
plot(30, error(30),'o','Color', TOLC(10,:),'LineWidth',1.5)
plot(60, error(60),'o','Color', TOLC(10,:),'LineWidth',1.5)
plot(305, error(305),'o','Color', TOLC(10,:),'LineWidth',1.5)
h_1=text(30,error(30),txt_1);
h_2=text(60,error(60),txt_2);
% h_3=text(305,error(305),txt_3);
 h_3=text(305,error(299),txt_3);
 h_4=text(300,error(280),txt_4);
set([h_1,h_2,h_3,h_4], 'Color', TOLC(10,:),'interpreter','latex','fontsize',11,'FontName', 'Times New Roman','fontweight','bold')
xlabel('R','fontsize',14,'interpreter','latex','FontName', 'Times New Roman','fontweight','bold')
ylabel('Error$$=\sqrt{\sum_{i=R+1}^{K} \sigma_{i}^{2}}$$','fontsize',14,'interpreter','latex','FontName', 'Times New Roman','fontweight','bold')
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
hold off



% r=30;
% Ur = U(:,1:r);
% Sr = S(1:r,1:r);
% Vr = V(:,1:r)'; % Truncate U,S,V using the rank r
% Xr = Ur*Sr*Vr; % Truncated matrix
% modeloutput_truncated=Xr;
% % % % save('SWE_POD_r20.mat','modeloutput_truncated','x','y','H','dt','t_save')
% save('SWE_POD_r30_new.mat','modeloutput_truncated','x','y','H','dt','t_save')
