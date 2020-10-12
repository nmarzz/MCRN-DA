clear;clc;

%Runs
%Q=1E-2, R=1E-2
epsQ=1E-1;
epsR=R;
Num=10;
Mult=40/Num;
for j=1:Num
numModes_physical=j*Mult;
[Time(:,j),RMSEsave(:,j), RMSEsave_proj(:,j), XCsave(:,j), XCprojsave(:,j), ESSsave(:,j), ResampPercent(:,j)]=PFL96run(numModes_physical,epsQ);
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

figure(1)
ObsErr = linspace(sqrt(epsR),sqrt(epsR),Num);
subplot(1,4,1)
plot(numModes,RMSEave,'b-')
hold on
plot(numModes,RMSEave_proj,'g-.')
plot(numModes,ObsErr,'r--')
legend('RMSE','RMSE(projected)','Obs Err')
xlabel('Model Dimension');
ylabel('RMSE');
title('RMSE');
hold off
subplot(1,4,2)
plot(numModes,XCave)
xlabel('Model Dimension');
ylabel('Mean Pattern Correlation Coefficient');
title('Mean Correlation Coeff');

subplot(1,4,3)
plot(numModes,ESSave)
xlabel('Model Dimension');
ylabel('Mean Effective Sample Size');
title('Mean Effective Sample Size');

subplot(1,4,4)
plot(numModes,RSpercent)
xlabel('Model Dimension');
ylabel('Resampling Percent');
title('Resampling Percent');

sgtitle(['Covariances: Q = ',num2str(epsQ), 'I, R = ', num2str(epsR), 'I']);

