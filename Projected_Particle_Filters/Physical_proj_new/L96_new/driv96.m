clear;clc;

%Runs
%Q=1E-2, R=1E-2
PhysicalProjection = 3;
DataProjection= 1;
inth = 1;
iOPPF = 1;


epsQ= 1;
epsR=1e-2;
Num=20;
Mult=40/Num;
for j=1:Num
    j
numModes_physical=j*Mult;
tolerance_physical = j*Mult;
[Time(:,j),RMSEsave(:,j), RMSEsave_proj(:,j), XCsave(:,j), XCprojsave(:,j), ESSsave(:,j), ResampPercent(:,j)]=PFL96run(tolerance_physical,numModes_physical,PhysicalProjection,DataProjection,epsQ,epsR);
end


filename = sprintf('PODp_%2d%2d_%4d_%2d_%4d.mat',PhysicalProjection,DataProjection,epsQ,epsR,inth)
params.PhysicalProjection = PhysicalProjection;
params.DataProjection = DataProjection;
params.epsQ = epsQ;
params.epsR=epsR;
params.iOPPF = iOPPF;
% params.ObsErr=ObsErr;
params.numModes=Mult:Mult:Num;
params.Num=Num;
params.Mult = Mult;
results.RMSEsave = RMSEsave;
results.RMSEsave_proj = RMSEsave_proj;
results.XCsave = XCsave; 
results.XCsave_proj =XCprojsave;

% % 
results.ESSsave= ESSsave;
results.ResampPercent =ResampPercent;
results.Time =Time;
save(filename,'params','results');
















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

% figure(1)
% ObsErr = linspace(sqrt(epsR),sqrt(epsR),Num);
% subplot(1,4,1)
% plot(numModes,RMSEave,'b-')
% hold on
% plot(numModes,RMSEave_proj,'g-.')
% plot(numModes,ObsErr,'r--')
% legend('RMSE','RMSE(projected)','Obs Err')
% xlabel('Model Dimension');
% ylabel('RMSE');
% title('RMSE');
% hold off
% subplot(1,4,2)
% plot(numModes,XCave)
% xlabel('Model Dimension');
% ylabel('Mean Pattern Correlation Coefficient');
% title('Mean Correlation Coeff');
% 
% subplot(1,4,3)
% plot(numModes,ESSave)
% xlabel('Model Dimension');
% ylabel('Mean Effective Sample Size');
% title('Mean Effective Sample Size');
% 
% subplot(1,4,4)
% plot(numModes,RSpercent)
% xlabel('Model Dimension');
% ylabel('Resampling Percent');
% title('Resampling Percent');
% 
% sgtitle(['Covariances: Q = ',num2str(epsQ), 'I, R = ', num2str(epsR), 'I']);
% 
