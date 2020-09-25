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

