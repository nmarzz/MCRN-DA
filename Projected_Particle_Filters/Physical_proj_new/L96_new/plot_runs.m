close all;clear all;clc;
%%
%  a_1=load('L968pod_ 1 1_1.000000e-01_1.000000e-02_   1.mat');%senario=1 (POD D; POD P )
% %a_1=load('L968dmd_ 2 1_   1_1.000000e-02_   1.mat');%senario=1 (DMD D; POD P)
% a_1=load('L968pod_ 1 1_   1_1.000000e-02_   1.mat');%
a_1=load('L963.5pod_ 1 1_   1_1.000000e-02_   1.mat');%senario=1 (POD D; POD P)
b_1=a_1.results;
c_1=a_1.params;
%  a_2=load('L968dmd_ 2 2_1.000000e-01_1.000000e-02_   1.mat');%senario=2 (DMD)
% a_2=load('L968dmd_ 2 2_   1_1.000000e-02_   1.mat');%senario=2 (DMD D, DMD P)
a_2=load( 'L963.5dmd_ 2 2_   1_1.000000e-02_   1.mat');%senario=2 (AUS D; DMD P)
% a_2=load('L968_ 2 1_   1_1.000000e-02_   1.mat');
b_2=a_2.results;
c_2=a_2.params;
% a_3=load('L968Aus_ 3 3_1.000000e-01_1.000000e-02_   1.mat');%senario=3 (AUS)
% a_3=load('L968Aus_ 3 3_   1_1.000000e-02_   1.mat');%senario=1 (DMD D; AUS P)
a_3=load('L963.5Aus_ 3 3_   1_1.000000e-02_   1.mat');%senario=2 (AUS D; AUS P)
% a_3= load('L968_ 3 1_   1_1.000000e-02_   1.mat');
b_3=a_3.results;
c_3=a_3.params;
%  a_0=load('L968_ 0 0_1.000000e-01_1.000000e-02_   1.mat');
a_0=load('L963.5_ 0 0_   1_1.000000e-02_   1.mat');
aa=a_0.results;
%%
PhysicalProjection=c_1.PhysicalProjection;
DataProjection=c_1.DataProjection;
epsQ=c_1.epsQ;
epsR=c_1.epsR;
% epsOmega_1=c_1.epsOmega;
iOPPF=c_1.iOPPF ;
numModes=c_1.numModes;
Num=c_1.Num;
%%
RMSEsave_1=b_1.RMSEsave;
RMSEsave_proj_1=b_1.RMSEsave_proj ;
XCsave_1=b_1.XCsave;
XCprojsave_1=b_1.XCsave_proj ;
ESSsave_1=b_1.ESSsave;
RSpercent_1=b_1.ResampPercent;
Time=b_1.Time;
% epsOmega=c_1.epsOmega;
%%
RMSEsave_2=b_2.RMSEsave;
RMSEsave_proj_2=b_2.RMSEsave_proj ;
XCsave_2=b_2.XCsave;
XCprojsave_2=b_2.XCsave_proj ;
ESSsave_2=b_2.ESSsave;
RSpercent_2=b_2.ResampPercent;
%%
RMSEsave_3=b_3.RMSEsave;
RMSEsave_proj_3=b_3.RMSEsave_proj ;
XCsave_3=b_3.XCsave;
XCprojsave_3=b_3.XCsave_proj ;
ESSsave_3=b_3.ESSsave;
RSpercent_3=b_3.ResampPercent;
%% 
RMSEsave_01=aa.RMSEsave;
RMSEsave_proj_01=aa.RMSEsave_proj ;
XCsave_01=aa.XCsave;
XCprojsave_01=aa.XCsave_proj ;
ESSsave_01=aa.ESSsave;
RSpercent_01=aa.ResampPercent;

%%
Mult = a_1.params.Mult;
for j=1:Num
    numModes(j)=j*Mult;
    RMSEave_1(j) = mean(RMSEsave_1(:,j));
    RMSEave_proj_1(j) = mean(RMSEsave_proj_1(:,j));
    XCave_1(j) = mean(XCsave_1(:,j));
    XCave_proj_1(j) = mean(XCprojsave_1(:,j));
    ESSave_1(j)=mean(ESSsave_1(:,j));
    RSpercent_1(j)=RSpercent_1(:,j);
end
%%
for j=1:Num
    numModes(j)=j*Mult;
    RMSEave_2(j) = mean(RMSEsave_2(:,j));
    RMSEave_proj_2(j) = mean(RMSEsave_proj_2(:,j));
    XCave_2(j) = mean(XCsave_2(:,j));
    XCave_proj_2(j) = mean(XCprojsave_2(:,j));
    ESSave_2(j)=mean(ESSsave_2(:,j));
    RSpercent_2(j)=RSpercent_2(:,j);
end
%%
for j=1:Num
    numModes(j)=j*Mult;
    RMSEave_3(j) = mean(RMSEsave_3(:,j));
    RMSEave_proj_3(j) = mean(RMSEsave_proj_3(:,j));
    XCave_3(j) = mean(XCsave_3(:,j));
    XCave_proj_3(j) = mean(XCprojsave_3(:,j));
    ESSave_3(j)=mean(ESSsave_3(:,j));
    RSpercent_3(j)=RSpercent_3(:,j);
end
%%
for j=1:Num
    numModes(j)=j*Mult;
    RMSEave_01(j) = mean(RMSEsave_01(:,j));
    RMSEave_proj_01(j) = mean(RMSEsave_proj_01(:,j));
    XCave_01(j) = mean(XCsave_01(:,j));
    XCave_proj_01(j) = mean(XCprojsave_01(:,j));
    ESSave_01(j)=mean(ESSsave_01(:,j));
    RSpercent_01(j)=RSpercent_01(:,j);
end
%%
TOLC=ptc12(12);
tt = tiledlayout(1,2);
ax(1) =nexttile;
plot(numModes,real(ESSave_01),'Color', TOLC(7,:),'LineStyle','-.','LineWidth', 1.5)
grid on
hold on
plot(numModes,real(ESSave_1),'Color', TOLC(1,:),'LineStyle','-','LineWidth', 1.5)
% grid on
% hold on
plot(numModes,real(ESSave_2),'Color', TOLC(11,:),'LineStyle',':','LineWidth', 1.5)
plot(numModes,real(ESSave_3),'Color', TOLC(4,:),'LineStyle','--','LineWidth', 1.5)
xlim([min(numModes) max(numModes)])
ylabel('Mean Effective Sample Size',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',12,...
    'FontName','Times')
title('Effective Sample Size',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',12,...
    'FontName','Times')
set(gca,...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',12,...
    'FontName','Times')

 hold off
ax(2)=nexttile;

plot(numModes,RSpercent_01,'Color', TOLC(7,:),'LineStyle','-.','LineWidth', 1.5)
grid on
hold on
plot(numModes,RSpercent_1,'Color', TOLC(1,:),'LineStyle','-','LineWidth', 1.5)
% grid on
% hold on
plot(numModes,RSpercent_2,'Color', TOLC(11,:),'LineStyle',':','LineWidth', 1.5)
plot(numModes,RSpercent_3,'Color', TOLC(4,:),'LineStyle','--','LineWidth', 1.5)
xlim([min(numModes) max(numModes)])

ylabel('Mean Resampling Percent',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',12,...
    'FontName','Times')
title('Resampling Percent',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',12,...
    'FontName','Times')
set(gca,...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',12,...
    'FontName','Times')
legend('Identity Projection','POD','DMD','AUS', ...
    'FontWeight','normal',...
    'FontSize',12,...
    'FontName','Times',...
    'Location','Best')
 hold off
sgtitle(['Covariances: $$Q = $$',num2str(epsQ), 'I, $$R = $$', num2str(epsR), 'I'],...
        'interpreter','latex',...
     'FontWeight','bold',...
    'FontSize',12,...
    'FontName','Times')

han=axes(gcf,'visible','off'); 
han.XLabel.Visible='on';
xlabel(han,'Model Dimensions',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',12,...
    'FontName','Times');
% print -depsc 8poddmdausQ-1R1e-2_essResamp%graph name

%%
figure(2)
TOLC=ptc12(12);
tt = tiledlayout(1,2);
ax(1) =nexttile;
plot(numModes,real(RMSEave_01),'Color', TOLC(7,:),'LineStyle','-.','LineWidth', 1.5)
grid on
hold on
plot(numModes,real(RMSEave_1),'Color', TOLC(1,:),'LineStyle','-','LineWidth', 1.5)
% grid on
% hold on
plot(numModes,real(RMSEave_2),'Color', TOLC(11,:),'LineStyle',':','LineWidth', 1.5)
plot(numModes,real(RMSEave_3),'Color', TOLC(4,:),'LineStyle','--','LineWidth', 1.5)
xlim([min(numModes) max(numModes)])
% xlabel('Model Dimension',...
%     'FontUnits','points',...
%     'FontWeight','normal',...
%     'FontSize',12,...
%     'FontName','Times')
ylabel('Mean RMSE',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',12,...
    'FontName','Times')
title('RMSE',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',12,...
    'FontName','Times')
set(gca,...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',12,...
    'FontName','Times')
% legend('Scenario(i)','Scenario(ii)','Scenario(iii)',... %scenario
%     'FontWeight','normal',...
%     'FontSize',12,...
%     'FontName','Times',...
%     'Location','Best')
% legend('boxoff')
 hold off
ax(2)=nexttile;
plot(numModes,real(XCave_01),'Color', TOLC(7,:),'LineStyle','-.','LineWidth', 1.5)
grid on
hold on
plot(numModes,real(XCave_1),'Color', TOLC(1,:),'LineStyle','-','LineWidth', 1.5)
% grid on
% hold on
plot(numModes,real(XCave_2),'Color', TOLC(11,:),'LineStyle',':','LineWidth', 1.5)
plot(numModes,real(XCave_3),'Color', TOLC(4,:),'LineStyle','--','LineWidth', 1.5)
xlim([min(numModes) max(numModes)])
% xlabel('Model Dimension',...
%     'FontUnits','points',...
%     'FontWeight','normal',...
%     'FontSize',12,...
%     'FontName','Times')
ylabel('Mean Pattern Correlation',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',12,...
    'FontName','Times')
title('Pattern Correlation',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',12,...
    'FontName','Times')
set(gca,...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',12,...
    'FontName','Times')
legend('Identity Projection','POD','DMD','AUS',...
    'FontWeight','normal',...
    'FontSize',12,...
    'FontName','Times',...
    'Location','Best')
 hold off
sgtitle(['Covariances: $$Q = $$',num2str(epsQ), 'I, $$R = $$', num2str(epsR), 'I'],...
        'interpreter','latex',...
     'FontWeight','bold',...
    'FontSize',12,...
    'FontName','Times')
% % , 'I, $\omega$ =', num2str(epsOmega)
% tt.TileSpacing = 'compact';
% tt.Padding = 'compact';
han=axes(gcf,'visible','off'); 
han.XLabel.Visible='on';
xlabel(han,'Model Dimensions',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',12,...
    'FontName','Times');
%% Here you can save the file as eps format which is good with Latex
% print -depsc SWE_physical_proj_DMD%graph name
% print -depsc SWE_physical_proj_DMD_inth1000%graph name
% print -depsc SWE_physical_proj_DMDPOD_inth1000%graph name
% print -depsc SWE_data_proj_POD_inth1%graph name
%print -depsc SWE_data_proj_dmd_inth1%graph name
%  print -depsc l968poddmdausQ1R1e-2_rmsexc%graph name
%  print -depsc l968poddmdausQ-1R1e-2_rmsexc%graph name
%  print -depsc l963poddmdausQ1R1e-2_rmsexc%graph name
%  print -depsc l963poddmdausQ-1R1e-2_rmsexc%graph name