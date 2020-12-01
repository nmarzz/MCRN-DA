close all; clear;clc;

%Runs
epsOmega =0.00001; %For inth = 1 %For inth = 1
%Observe every inth variable.
inth=1;
Numsteps=100;
%Q=1E-2, R=1E-2
epsQ=1;
epsR=1E-2;
Num=10;
Mult=100/Num;

%% Type of particle filter
% Use of standard PF or OP-PF (iOPPF=0 => standard PF, iOPPF=1 => OP-PF)
iOPPF=1; 
%% Projection_type(0 = no projection, 1 POD, 2 DMD, 3 AUS)
PhysicalProjection =1;
DataProjection =1  ;
% Observed variables scenario following Paulina et al.
% scenario 1: observe inth u and v's
% scenario 2: observe inth everything
% scenario 3: observe inth h
scenario =3 ;
for j=1:Num
numModes_physical=j*Mult;%DMD
tolerance_physical=j*Mult;%POD
[Time(:,j),RMSEsave(:,j), RMSEsave_proj(:,j), XCsave(:,j), XCprojsave(:,j), ESSsave(:,j), ResampPercent(:,j)]=PFSWErun...
    (numModes_physical,epsQ,epsR,tolerance_physical,iOPPF,PhysicalProjection,DataProjection,scenario,epsOmega, inth,Numsteps);
end

% %For Plots
% for j=1:Num
% numModes(j)=j*Mult;
% RMSEave(j) = mean(RMSEsave(:,j));
% RMSEave_proj(j) = mean(RMSEsave_proj(:,j));
% XCave(j) = mean(XCsave(:,j));
% XCave_proj(j) = mean(XCprojsave(:,j));
% ESSave(j)=mean(ESSsave(:,j));
% RSpercent(j)=ResampPercent(:,j);
% end
% TOLC=ptc12(9);
% tt = tiledlayout(1,4);
% ObsErr = linspace(sqrt(epsR),sqrt(epsR),Num);
% nexttile
% plot(numModes,RMSEave,'Color', TOLC(1,:),'LineStyle','--','LineWidth', 2)
% hold on
% plot(numModes,RMSEave_proj,'Color', TOLC(7,:),'LineStyle','--','Marker','+','LineWidth', 2)
% % semilogy(numModes,RMSE_no_proj,'Color', TOLC(7,:),'LineStyle',':','Marker','.','LineWidth',2)
% plot(numModes,ObsErr,'k-.','LineWidth', 2)
% grid on
% legend('Model Space','Projected Space','Observation Error','Location', 'Best','fontsize',11,'interpreter','latex','FontName', 'Times New Roman','fontweight','bold')
% xlabel('Model Dimension','fontsize',14,'interpreter','latex','FontName', 'Times New Roman','fontweight','bold')
% ylabel('Mean RMSE','fontsize',14,'interpreter','latex','FontName', 'Times New Roman','fontweight','bold')
% title('RMSE','FontName', 'Times New Roman', 'FontSize', 14)
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
% hold off
% 
% nexttile
% plot(numModes,real(XCave),'Color', TOLC(1,:),'LineStyle','--','LineWidth', 2)
% hold on
% plot(numModes, real(XCave_proj),'Color', TOLC(7,:),'LineStyle',':','Marker','.','LineWidth',2)
% grid on
% xlabel('Model Dimension','fontsize',14,'interpreter','latex','FontName', 'Times New Roman','fontweight','bold')
% ylabel('Mean Pattern Correlations ','fontsize',14,'interpreter','latex','FontName', 'Times New Roman','fontweight','bold')
% legend('Model Space','Projected Space','Location', 'Best','fontsize',11,'interpreter','latex','FontName', 'Times New Roman','fontweight','bold')
% title('Mean Correlation Coeff','FontName', 'Times New Roman', 'FontSize', 14);
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
% hold off
% 
% nexttile
% plot(numModes,ESSave,'Color', TOLC(4,:),'LineStyle','-')
% grid on
% xlabel('Model Dimension','fontsize',14,'interpreter','latex','FontName', 'Times New Roman','fontweight','bold');
% ylabel('Mean Effective Sample Size','fontsize',14,'interpreter','latex','FontName', 'Times New Roman','fontweight','bold');
% title('Mean Effective Sample Size','FontName', 'Times New Roman', 'FontSize', 14);
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
% 
% nexttile
% plot(numModes,RSpercent,'Color', TOLC(4,:),'LineStyle','-')
% grid on
% xlabel('Model Dimension','fontsize',14,'interpreter','latex','FontName', 'Times New Roman','fontweight','bold');
% ylabel('Resampling Percent','fontsize',14,'interpreter','latex','FontName', 'Times New Roman','fontweight','bold');
% title('Resampling Percent','FontName', 'Times New Roman', 'FontSize', 14);
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
% 
% sgtitle(['Covariances: Q = ',num2str(epsQ), 'I, R = ', num2str(epsR), 'I'],'fontsize',14,'interpreter','latex','FontName', 'Times New Roman','fontweight','bold');
% tt.TileSpacing = 'compact';
% tt.Padding = 'compact';
% %%
% figure(2)
% plot(Time,real(ESSsave))
% grid on
% xlabel('Time','fontsize',12,'interpreter','latex','FontName', 'Times New Roman','fontweight','bold');
% ylabel('ESS','fontsize',12,'interpreter','latex','FontName', 'Times New Roman','fontweight','bold');
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)

L=20;
% % %% Save to mat file
filename = sprintf('SWEp_%2d%2d_%2d_%2d_%4d_%d_%d.mat',PhysicalProjection,DataProjection,epsQ,epsR,inth,scenario,L)
params.PhysicalProjection = PhysicalProjection;
params.DataProjection = DataProjection;
params.epsQ = epsQ;
params.epsR=epsR;
params.iOPPF = iOPPF;
% params.ObsErr=ObsErr;
params.numModes=numModes_physical;
params.Num=Num;
params.epsOmega=epsOmega;

results.RMSEsave = RMSEsave;
results.RMSEsave_proj = RMSEsave_proj;

results.XCsave = XCsave; 
results.XCsave_proj =XCprojsave;
% % 
results.ESSsave= ESSsave;

results.ResampPercent =ResampPercent;
results.Time =Time;
save(filename,'params','results');
