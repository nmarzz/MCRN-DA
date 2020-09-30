clear;clc;

%Runs
%Q=1E-2, R=1E-2
epsQ=1;
epsR=1E-2;
Num=10;
Mult=100/Num;

%% Type of particle filter
% Use of standard PF or OP-PF (iOPPF=0 => standard PF, iOPPF=1 => OP-PF)
iOPPF=1;
%% Projection_type(0 = no projection, 1 POD, 2 DMD, 3 AUS)
PhysicalProjection =2;
DataProjection =1;
for j=1:Num
numModes_physical=j*Mult;%DMD
tolerance_physical=j*Mult;%POD
[Time(:,j),RMSEsave(:,j), RMSEsave_proj(:,j), XCsave(:,j), XCprojsave(:,j), ESSsave(:,j), ResampPercent(:,j)]=PFSWErun...
    (numModes_physical,epsQ,tolerance_physical,iOPPF,PhysicalProjection,DataProjection);
end

%For Plots
for j=1:Num
numModes(j)=j*Mult;
RMSEave(j) = mean(RMSEsave(:,j));
RMSEave_proj(j) = mean(RMSEsave_proj(:,j));
XCave(j) = mean(XCsave(:,j));
XCave_proj(j) = mean(XCprojsave(:,j));
ESSave(j)=mean(ESSsave(:,j));
RSpercent(j)=ResampPercent(:,j);
end
% TOLC=ptc12(9);
% tt = tiledlayout(1,4);
% ObsErr = linspace(sqrt(epsR),sqrt(epsR),Num);
% nexttile
% plot(numModes,RMSEave,'Color', TOLC(1,:),'LineStyle','-','LineWidth', 2)
% hold on
% plot(numModes,RMSEave_proj,'Color', TOLC(7,:),'LineStyle','--','Marker','+','LineWidth', 2)
% % semilogy(numModes,RMSE_no_proj,'Color', TOLC(7,:),'LineStyle',':','Marker','.','LineWidth',2)
% plot(numModes,ObsErr,'k-.','LineWidth', 2)
% grid on
% legend('Model Space','Projected Space','Observation Error','Location', 'Best','fontsize',10,'interpreter','latex','FontName', 'Times New Roman','fontweight','bold')
% xlabel('Model Dimension','fontsize',12,'interpreter','latex','FontName', 'Times New Roman','fontweight','bold')
% ylabel('Mean RMSE','fontsize',12,'interpreter','latex','FontName', 'Times New Roman','fontweight','bold')
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)
% hold off
% 
% nexttile
% plot(numModes,real(XCave),'Color', TOLC(1,:),'LineStyle','-','LineWidth', 2)
% hold on
% plot(numModes, real(XCave_proj),'Color', TOLC(7,:),'LineStyle',':','Marker','.','LineWidth',2)
% grid on
% xlabel('Model Dimension','fontsize',12,'interpreter','latex','FontName', 'Times New Roman','fontweight','bold')
% ylabel('Mean Pattern Correlations ','fontsize',12,'interpreter','latex','FontName', 'Times New Roman','fontweight','bold')
% legend('Model Space','Projected Space','Location', 'Best','fontsize',10,'interpreter','latex','FontName', 'Times New Roman','fontweight','bold')
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)
% hold off
% 
% nexttile
% plot(numModes,ESSave)
% grid on
% xlabel('Model Dimension','fontsize',12,'interpreter','latex','FontName', 'Times New Roman','fontweight','bold');
% ylabel('Mean Effective Sample Size','fontsize',12,'interpreter','latex','FontName', 'Times New Roman','fontweight','bold');
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)
% 
% nexttile
% plot(numModes,RSpercent)
% grid on
% xlabel('Model Dimension','fontsize',12,'interpreter','latex','FontName', 'Times New Roman','fontweight','bold');
% ylabel('Resampling Percent','fontsize',12,'interpreter','latex','FontName', 'Times New Roman','fontweight','bold');
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)
% 
% sgtitle(['Covariances: Q = ',num2str(epsQ), 'I, R = ', num2str(epsR), 'I'],'fontsize',14,'interpreter','latex','FontName', 'Times New Roman','fontweight','bold');
% tt.TileSpacing = 'compact';
% tt.Padding = 'compact';
% %%
% figure(2)
% plot(Time,ESSsave)
% grid on
% xlabel('Time','fontsize',12,'interpreter','latex','FontName', 'Times New Roman','fontweight','bold');
% ylabel('ESS','fontsize',12,'interpreter','latex','FontName', 'Times New Roman','fontweight','bold');
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)
%% Save to mat file
filename = sprintf('SWE_%2d%2d_%4d_%2d.mat',PhysicalProjection,DataProjection,epsQ,epsR)
params.PhysicalProjection = PhysicalProjection;
params.DataProjection = DataProjection;
params.epsQ = epsQ;
params.epsR=epsR;
params.iOPPF = iOPPF;
params.ObsErr=ObsErr;
params.numModes_physical=numModes_physical;

results.RMSEave_orig = RMSEave;
results.RMSEave_proj = RMSEave_proj;
results.RMSEsave = RMSEsave;
results.RMSEsave_proj = RMSEsave_proj;

results.XCave= XCave;
results.XCave_proj = XCave_proj;
results.XCsave = XCsave;
results.XCsave_proj =XCprojsave;

results.ESSave= ESSave;
results.ESSsave= ESSsave;

results.ResampPercent =RSpercent;
results.Time =Time;
save(filename,'params','results');