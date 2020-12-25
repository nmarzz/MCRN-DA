close all;clear all;clc;
%% N=1000
% a_1=load('L968p_ 1 1_1.000000e-02_1.000000e-02_1000.mat');
% a_2=load('L968p_ 1 1_1.000000e-01_1.000000e-02_1000.mat');

%% Q=1, POD, N=400,F=3,  F=4, F=6, F=8
% a_1=load('L963p_ 1 1_ 1_1.000000e-02_ 400.mat');
% a_2= load('L964p_ 1 1_ 1_1.000000e-02_ 400.mat');
% a_3= load('L966p_ 1 1_ 1_1.000000e-02_ 400.mat');
% a_4= load('L968p_ 1 1_ 1_1.000000e-02_ 400.mat');

%% Q=1E-1, POD, N=400, F=3,F=4, F=6, F=8
% a_1=load('L963p_ 1 1_1.000000e-01_1.000000e-02_ 400.mat');%F=3
% a_2=load ('L964p_ 1 1_1.000000e-01_1.000000e-02_ 400.mat');%F=4
% a_3=load('L966p_ 1 1_1.000000e-01_1.000000e-02_ 400.mat');%F=6
% a_4=load('L968p_ 1 1_1.000000e-01_1.000000e-02_ 400.mat');%F=8

%% POD and DMD Q=1, F=6,8, N=400
a_1=load('L966p_ 1 1_ 1_1.000000e-02_ 400.mat');% POD, F=6
a_2=load('L966p_ 2 2_ 1_1.000000e-02_ 400.mat'); %DMD, F=6
a_3= load('L968p_ 1 1_ 1_1.000000e-02_ 400.mat');%POD, F=8
a_4=load('L968p_ 2 2_ 1_1.000000e-02_ 400.mat'); %DMD, F=8

%% POD and DMD Q=1E-1, F=6,8, N=400
% a_1=load('L966p_ 1 1_1.000000e-01_1.000000e-02_ 400.mat');% POD, F=6
% a_2=load('L966p_ 2 2_1.000000e-01_1.000000e-02_ 400.mat'); %DMD, F=6
% a_3=load('L968p_ 1 1_1.000000e-01_1.000000e-02_ 400.mat');%POD, F=8
% a_4=load('L968p_ 2 2_1.000000e-01_1.000000e-02_ 400.mat'); %DMD, F=8


%% POD and DMD Q=1E-1, F=3,4, N=400
% a_1=load('L963p_ 1 1_1.000000e-01_1.000000e-02_ 400.mat');%POD, F=3
% a_2= load('L963p_ 2 2_1.000000e-01_1.000000e-02_ 400.mat');%DMD,F=3
% a_3= load('L964p_ 1 1_1.000000e-01_1.000000e-02_ 400.mat');%POD, F=4
% a_4=  load('L964p_ 2 2_1.000000e-01_1.000000e-02_ 400.mat');%DMD,F=4

%% POD and DMD Q=1, F=3,4, N=400
% a_1= load('L963p_ 1 1_ 1_1.000000e-02_ 400.mat');%POD, F=3
% a_2= load('L963p_ 2 2_ 1_1.000000e-02_ 400.mat');%DMD,F=3
% a_3=load('L964p_ 1 1_ 1_1.000000e-02_ 400.mat');%POD, F=4
% a_4=  load('L964p_ 2 2_ 1_1.000000e-02_ 400.mat');%DMD,F=4

%%
b_1=a_1.results;
c_1=a_1.params;
%%
b_2=a_2.results;
c_2=a_2.params;
%%
b_3=a_3.results;
c_3=a_3.params;
%%
b_4=a_4.results;
c_4=a_4.params;
%%
PhysicalProjection=c_1.PhysicalProjection;
DataProjection=c_1.DataProjection;
epsQ=c_1.epsQ;
epsR=c_1.epsR;
numModes=c_1.numModes;
Num=c_1.Num;
Mult=c_1.Mult;
%%
RMSEsave_1=b_1.RMSEsave;
RSpercent_1=b_1.ResampPercent;
Time=b_1.Time;
%%
RMSEsave_2=b_2.RMSEsave;
RSpercent_2=b_2.ResampPercent;
%%
RMSEsave_3=b_3.RMSEsave;
RSpercent_3=b_3.ResampPercent;
%%
RMSEsave_4=b_4.RMSEsave;
RSpercent_4=b_4.ResampPercent;
%%
for j=1:Num
    numModes(j)=j*Mult;
    RMSEave_1(j) = mean(RMSEsave_1(:,j));
    RSpercent_1(j)=RSpercent_1(:,j);
    RMSEave_2(j) = mean(RMSEsave_2(:,j));
    RSpercent_2(j)=RSpercent_2(:,j);
    RMSEave_3(j) = mean(RMSEsave_3(:,j));
    RSpercent_3(j)=RSpercent_3(:,j);
    RMSEave_4(j) = mean(RMSEsave_4(:,j));
    RSpercent_4(j)=RSpercent_4(:,j);
end
%%
TOLC=ptc12(12);
figure(1)
tt = tiledlayout(1,2);
ObsErr = linspace(sqrt(epsR),sqrt(epsR),Num);
nexttile
plot(numModes,real(RMSEave_1),'Color', TOLC(1,:),'LineStyle','-','LineWidth', 1.5)
grid on
hold on
plot(numModes,real(RMSEave_2),'Color', TOLC(7,:),'LineStyle',':','LineWidth', 1.5)
plot(numModes,real(RMSEave_3),'Color', TOLC(4,:),'LineStyle','--','LineWidth', 1.5)
plot(numModes,real(RMSEave_4),'Color', TOLC(11,:),'LineStyle','-.','LineWidth', 1.5)
plot(numModes,ObsErr,'k-.','LineWidth', 1.5)
% legend('$Q = 1e-2$','$Q = 1e-1$',...
%     'interpreter','latex',...
%     'FontWeight','normal',...
%     'FontSize',12,...
%     'FontName','Times',...
%     'Location','Best')
%
% legend('F=3, POD','F=3, DMD','F=4, POD','F=4, DMD','ObsErr',...
%     'FontWeight','normal',...
%     'FontSize',12,...
%     'FontName','Times',...
%     'Location','Best')
legend('F=6, POD','F=6, DMD','F=8, POD','F=8, DMD','ObsErr',...
    'FontWeight','normal',...
    'FontSize',12,...
    'FontName','Times',...
    'Location','Best')
% legend('POD','DMD','ObsErr',...
%     'FontWeight','normal',...
%     'FontSize',12,...
%     'FontName','Times',...
%     'Location','Best')
xlabel('Model Dimensions',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',14,...
    'FontName','Times')
ylabel('Mean RMSE',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',14,...
    'FontName','Times')

title('RMSE',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',14,...
    'FontName','Times')
set(gca,...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',14,...
    'FontName','Times')
xlim([min(numModes) max(numModes)])
set(gca,...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',14,...
    'FontName','Times')
hold off
%%
nexttile
plot(numModes,RSpercent_1,'Color', TOLC(1,:),'LineStyle','-','LineWidth', 1.5)
grid on
hold on
plot(numModes,RSpercent_2,'Color', TOLC(7,:),'LineStyle',':','LineWidth', 1.5)
plot(numModes,RSpercent_3,'Color', TOLC(4,:),'LineStyle','--','LineWidth', 1.5)
plot(numModes,RSpercent_4,'Color', TOLC(11,:),'LineStyle','-.','LineWidth', 1.5)
xlim([min(numModes) max(numModes)])
% legend('L=05','L=10','L=20',...
%     'FontWeight','normal',...
%     'FontSize',12,...
%     'FontName','Times',...
%     'Location','Best')
% legend('$Q = 1e-2$','$Q = 1e-1$',...
%     'interpreter','latex',...
%     'FontWeight','normal',...
%     'FontSize',12,...
%     'FontName','Times',...
%     'Location','Best')
% legend('F=3, POD','F=3, DMD','F=4, POD','F=4, DMD',...
%     'FontWeight','normal',...
%     'FontSize',12,...
%     'FontName','Times',...
%     'Location','Best')
legend('F=6, POD','F=6, DMD','F=8, POD','F=8, DMD',...
    'FontWeight','normal',...
    'FontSize',12,...
    'FontName','Times',...
    'Location','Best')
% legend('POD','DMD',...
%     'FontWeight','normal',...
%     'FontSize',12,...
%     'FontName','Times',...
%     'Location','Best')
xlabel('Model Dimensions',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',14,...
    'FontName','Times')

ylabel('Mean Resampling Percent',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',14,...
    'FontName','Times')
title('Resampling Percent',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',14,...
    'FontName','Times')
set(gca,...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',14,...
    'FontName','Times')
hold off
sgtitle(['Covariances: $Q = $',num2str(epsQ), 'I, $R = $', num2str(epsR), 'I'],...
    'interpreter','latex',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Times')
tt.TileSpacing = 'compact';
tt.Padding = 'compact';

%% Here you can save the file as eps format which is good with Latex
% print -depsc l968podQ1R1e-2%graph name