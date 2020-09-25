clear all;close all;clc
% load('SWE.mat')
load('SWE_run_4days.mat');
DataMatrix= x_save;%snapshot matrix
r=130;
VALUE=true;
M = mean(DataMatrix, 2); % 2 = mean across columns
stepVALUE=[1 10 30 60];
out1 = dmd( DataMatrix, dt, r, 'removemean', VALUE,'step',stepVALUE,'sortbyb', VALUE);
bM = norm(M);
Mnormalized = M/bM;
% manually adding the mean value of data back among the modes
out1.Phi = [Mnormalized(:), out1.Phi];
out1.b = [bM; out1.b];
out1.omega = [0; out1.omega];
out1.lambda = [1; out1.lambda];
%%
Omega=out1.omega ;
A=real(Omega)>0;
B=real(Omega)<0;
C=real(Omega)==0;
scatter(A,imag(Omega),'g','filled')
hold on; grid on;
scatter(B,imag(Omega),'b','filled')
plot(C(1,1),imag(Omega),'r','filled')
