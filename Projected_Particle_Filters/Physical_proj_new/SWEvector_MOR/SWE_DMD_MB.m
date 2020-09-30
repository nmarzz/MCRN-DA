clear all;close all;clc
load('SWE_run_4days.mat');
DataMatrix= x_save;%snapshot matrix
r=130;
VALUE=true;
M = mean(DataMatrix, 2); % 2 = mean across columns
stepVALUE=[1 10 30 60];
out = dmd( DataMatrix, dt, r, 'removemean', VALUE,'step',stepVALUE,'sortbyb', VALUE);
bM = norm(M);
Mnormalized = M/bM;
% manually adding the mean value of data back among the modes
out.Phi = [Mnormalized(:), out.Phi];
out.b = [bM; out.b];
out.omega = [0; out.omega];
out.lambda = [1; out.lambda];
%%
Omega=out.omega ;
unstable_index=real(Omega)>1E-6;%~not stable
stable_index=real(Omega)<-1E-6;

CC=abs(real(Omega))<1E-6;
figure(2)
scatter(real(Omega(unstable_index)),imag(Omega(unstable_index)),'g','filled')
hold on; grid on;
scatter(real(Omega(stable_index)),imag(Omega(stable_index)),'b','filled')
scatter(real(Omega(CC)),imag(Omega(CC)),'r','filled')
%%
Lambda=out.lambda;
figure(3)
scatter(real(Lambda), imag(out.lambda),'filled')
% figure(3)
% scatter(real(Lambda), imag(out1.lambda),1+50*abs(out1.b(:,1))/max(abs(out1.b(:,1))),'filled')
%%
bb=out.b;
Phi=out.Phi(:,stable_index);
figure(4)
stem(imag(Lambda(2:end)), abs(out.meanL2norm),'filled')
 set(gca, 'yscale','log')