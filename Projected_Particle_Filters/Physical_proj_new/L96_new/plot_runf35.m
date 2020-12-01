close all;clear all;clc;
% %% F=3.5, Q=0.01, N=40
% a_1=load('L93.5_ 1 1_0.010000_1.000000e-02_ 1_40.mat');%POD 
% a_2=load('L93.5_ 2 2_0.010000_1.000000e-02_ 1_40.mat') ;%DMD
% a_3=load('L93.5_ 3 3_0.010000_1.000000e-02_ 1_40.mat'); %AUS
% %% F=3.5, Q=0.05, N=40
% a_1=load('L963.5_ 1 1_0.050000_1.000000e-02_ 1_40.mat');%POD 
% a_2=load('L963.5_ 2 2_0.050000_1.000000e-02_ 1_40.mat') ;%DMD
% a_3=load('L963.5_ 3 3_0.050000_1.000000e-02_ 1_40.mat'); %AUS
% %%  F=3.5, Q=0.1, N=40 
% a_1=load( 'L93.5_ 1 1_0.100000_1.000000e-02_ 1_40.mat');%POD 
% a_2=load('L93.5_ 2 2_0.100000_1.000000e-02_ 1_40.mat') ;%DMD
% a_3=load('L963.5_ 3 3_5.000000e-01_1.000000e-02_   1_40.mat'); %AUS
% %%  F=3.5, Q=0.5, N=40 
% a_1=load( 'L963.5_ 1 1_5.000000e-01_1.000000e-02_   1_40.mat');%POD 
% a_2=load('L963.5_ 2 2_5.000000e-01_1.000000e-02_   1_40.mat') ;%DMD
% a_3=load('L963.5_ 3 3_5.000000e-01_1.000000e-02_   1_40.mat'); %AUS
% %% F=3.5, Q=1, N=40
% a_1=load('L963.5_ 1 1_1.000000_1.000000e-02_ 1_40.mat');%POD 
% a_2=load('L963.5_ 2 2_1.000000_1.000000e-02_ 1_40.mat') ;%DMD
% a_3=load('L963.5_ 3 3_1.000000_1.000000e-02_ 1_40.mat'); %AUS
% %% F=3.5, Q=0.01, N=100
% a_1= load('L93.5_ 1 1_0.010000_1.000000e-02_ 1_100.mat');%POD 
% a_2= load('L93.5_ 2 2_0.010000_1.000000e-02_ 1_100.mat') ;%DMD
% a_3= load('L93.5_ 3 3_0.010000_1.000000e-02_ 1_100.mat'); %AUS
% %% F=3.5, Q=0.05, N=100
% a_1=load('L963.5_ 1 1_0.050000_1.000000e-02_ 1_100.mat');%POD 
% a_2=load('L963.5_ 2 2_0.050000_1.000000e-02_ 1_100.mat') ;%DMD
% a_3=load('L963.5_ 3 3_0.050000_1.000000e-02_ 1_100.mat'); %AUS
% %%  F=3.5, Q=0.1, N=100
% a_1=load('L93.5_ 1 1_0.100000_1.000000e-02_ 1_100.mat');%POD 
% a_2=load('L93.5_ 2 2_0.100000_1.000000e-02_ 1_100.mat');%DMD
% a_3=load('L93.5_ 3 3_0.100000_1.000000e-02_ 1_100.mat'); %AUS
% %%  F=3.5, Q=0.5, N=100
% a_1=load( 'L963.5_ 1 1_5.000000e-01_1.000000e-02_   1_100.mat');%POD 
% a_2=load('L963.5_ 2 2_5.000000e-01_1.000000e-02_   1_100.mat') ;%DMD
% a_3=load('L963.5_ 3 3_5.000000e-01_1.000000e-02_   1_100.mat'); %AUS
% %% F=3.5, Q=1, N=100 
% a_1=load('L963.5_ 1 1_1.000000_1.000000e-02_ 1_100.mat');%POD 
% a_2=load('L963.5_ 2 2_1.000000_1.000000e-02_ 1_100.mat') ;%DMD
% a_3=load('L963.5_ 3 3_1.000000_1.000000e-02_ 1_100.mat'); %AUS
%% iOPPF=0%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% F=3.5, Q=0.01, N=40
a_1=load('L93.5_ 1 1_0.010000_1.000000e-02_ 0_40.mat');%POD 
a_2=load('L93.5_ 2 2_0.010000_1.000000e-02_ 0_40.mat') ;%DMD
a_3=load('L93.5_ 3 3_0.010000_1.000000e-02_ 0_40.mat'); %AUS
% %% F=3.5, Q=1, N=40
% a_1=load( 'L93.5_ 1 1_1.000000_1.000000e-02_ 0_40.mat');%POD 
% a_2=load( 'L93.5_ 2 2_1.000000_1.000000e-02_ 0_40.mat') ;%DMD
% a_3=load('L93.5_ 3 3_1.000000_1.000000e-02_ 0_40.mat'); %AUS

%%
b_1=a_1.results;
c_1=a_1.params;

b_2=a_2.results;
c_2=a_2.params;

b_3=a_3.results;
c_3=a_3.params;
%%
PhysicalProjection=c_1.PhysicalProjection;
DataProjection=c_1.DataProjection;
epsQ=c_1.epsQ;
epsR=c_1.epsR;
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
TOLC=ptc12(12);
figure(2)
tt = tiledlayout(2,2);
ObsErr = linspace(sqrt(epsR),sqrt(epsR),Num);
nexttile
plot(numModes,real(RMSEave_1),'Color', TOLC(1,:),'LineStyle','-','LineWidth', 1.5)
grid on
hold on
plot(numModes,real(RMSEave_2),'Color', TOLC(11,:),'LineStyle',':','LineWidth', 1.5)
plot(numModes,real(RMSEave_3),'Color', TOLC(4,:),'LineStyle','--','LineWidth', 1.5)
plot(numModes,ObsErr,'k-.','LineWidth', 1.5)

legend('POD','DMD','AUS','ObsErr',...
    'FontWeight','normal',...
    'FontSize',12,...
    'FontName','Times',...
    'Location','Best')

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
hold off
%% 
nexttile
plot(numModes,real(ESSave_1),'Color', TOLC(1,:),'LineStyle','-','LineWidth', 1.5)
grid on
hold on
plot(numModes,real(ESSave_2),'Color', TOLC(11,:),'LineStyle',':','LineWidth', 1.5)
plot(numModes,real(ESSave_3),'Color', TOLC(4,:),'LineStyle','--','LineWidth', 1.5)
legend('POD','DMD','AUS',...
    'FontWeight','normal',...
    'FontSize',12,...
    'FontName','Times',...
    'Location','Best')

xlabel('Model Dimensions',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',14,...
    'FontName','Times')
ylabel('Mean Effective Sample Size',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',14,...
    'FontName','Times')

title('Effective Sample Size',...
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

%%
% nexttile([2 1])
 nexttile
plot(numModes,real(XCave_1),'Color', TOLC(1,:),'LineStyle','-','LineWidth', 1.5)
grid on
hold on
plot(numModes,real(XCave_2),'Color', TOLC(11,:),'LineStyle',':','LineWidth', 1.5)
plot(numModes,real(XCave_3),'Color', TOLC(4,:),'LineStyle','--','LineWidth', 1.5)
legend('POD','DMD','AUS',...
    'FontWeight','normal',...
    'FontSize',12,...
    'FontName','Times',...
    'Location','Best')

xlabel('Model Dimensions',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',14,...
    'FontName','Times')
ylabel('Mean Pattern Correlations ',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',14,...
    'FontName','Times')

title('Pattern Correlations',...
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
hold off
%%
nexttile

plot(numModes,RSpercent_1,'Color', TOLC(1,:),'LineStyle','-','LineWidth', 1.5)
grid on
hold on
plot(numModes,RSpercent_2,'Color', TOLC(11,:),'LineStyle',':','LineWidth', 1.5)
plot(numModes,RSpercent_3,'Color', TOLC(4,:),'LineStyle','--','LineWidth', 1.5)
xlim([min(numModes) max(numModes)])

legend('POD','DMD','AUS',...
    'FontWeight','normal',...
    'FontSize',12,...
    'FontName','Times',...
    'Location','Best')
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

% print -depsc l963poddmdausQ1R1e-2N40%graph name
% print -depsc l968poddmdausQ1R1e-2%graph name
% print -depsc l968poddmdausQ1e-1R1e-2%graph name
% print -depsc l968podQ1R1e-2%graph name