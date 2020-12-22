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
TOLC=ptc12(12);
figure(2)
semilogy(sig(1:350),'-ok','LineWidth',0.9), grid on, hold on
txt_1 =(  ' $\longrightarrow  r =20$');
txt_2 =( ' $\longrightarrow  r = 30$');
txt_3 = (' $\longrightarrow r=60$' );
txt_4 =(  ' $\nearrow $');
txt_5 =(  ' $ r= 305$');
plot(20, sig(20),'o','Color', TOLC(10,:),'LineWidth',1.5)
plot(30, sig(30),'o','Color', TOLC(10,:),'LineWidth',1.5)
plot(60, sig(60),'o','Color', TOLC(10,:),'LineWidth',1.5)
plot(305, sig(305),'o','Color', TOLC(10,:),'LineWidth',1.5)
% plot(130, sig(130),'o','Color', TOLC(10,:),'LineWidth',1.5)
h=text(20,sig(20),txt_1);
h_1=text(30,sig(30),txt_2);
h_3=text(60,sig(60),txt_3);
h_4=text(305,sig(290),txt_4);
h_5=text(302,sig(240),txt_5);
set([h,h_1,h_3,h_4,h_5], 'Color', TOLC(10,:),'interpreter','latex','fontsize',13,'FontName', 'Times','fontweight','bold')
xlim([0,350]);
xlabel('i',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',14,...
    'FontName','Times')
ylabel('Singular Values $\sigma_{i}^{2}$',...
    'interpreter','latex',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Times')
set(gca,...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',14,...
    'FontName','Times')
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
semilogy(error,'k-o','LineWidth',0.9),grid on,hold on;
txt_1 =(  ' $\longrightarrow  r =20$');
txt_2 =( ' $\longrightarrow  r = 30$');
txt_3 = (' $\longrightarrow r=60$' );
txt_4 =(  ' $\nearrow $');
txt_5 =(  ' $ r= 305$');
plot(20, error(20),'o','Color', TOLC(10,:),'LineWidth',1.5)
plot(30, error(30),'o','Color', TOLC(10,:),'LineWidth',1.5)
plot(60, error(60),'o','Color', TOLC(10,:),'LineWidth',1.5)
plot(305, error(305),'o','Color', TOLC(10,:),'LineWidth',1.5)
h=text(20,error(20),txt_1);
h_1=text(30,error(30),txt_2);
h_2=text(60,error(60),txt_3);
h_3=text(305,error(299),txt_4);
h_4=text(302,error(280),txt_5);
set([h,h_1,h_2,h_3,h_4], 'Color', TOLC(10,:),'interpreter','latex','fontsize',13,'FontName', 'Times','fontweight','bold')
% ylabel('Error$$=\sqrt{\sum_{i=R+1}^{K} \sigma_{i}^{2}}$$','fontsize',14,'interpreter','latex','FontName', 'Times New Roman','fontweight','bold')

xlabel('r',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',14,...
    'FontName','Times')
ylabel('Error',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',14,...
    'FontName','Times')
set(gca,...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',14,...
    'FontName','Times')

xlim([0,350]);
hold off

% r=30;
% Ur = U(:,1:r);
% Sr = S(1:r,1:r);
% Vr = V(:,1:r)'; % Truncate U,S,V using the rank r
% Xr = Ur*Sr*Vr; % Truncated matrix
% modeloutput_truncated=Xr;
% % % % save('SWE_POD_r20.mat','modeloutput_truncated','x','y','H','dt','t_save')
% save('SWE_POD_r30_new.mat','modeloutput_truncated','x','y','H','dt','t_save')
