clear all;close all;clc
load('SWE.mat')
% whos
modeloutput=x_save(:,1440:end);%snapshot matrix
t_save=t_save(:,1440:end);
L=length(t_save);%number of snapshots
m=length(y);%x-dimension
n=length(x);%y-dimension
N_gridpoints=m*n;%each snapshot has x*y-dimensions

%% DMD method
X1 =modeloutput(:,1:end-1);%The data matrix
X2 = modeloutput(:,2:end); %shifted data matrix

%% Step1:SVD to the data matrix
[U1,S1,V1] = svd(X1, 'econ');
r=10 ;%number of singular values that we will use
Ur=U1(:,1:r);Sr=S1(1:r,1:r);Vr=V1(:,1:r);%trancated matrices
Atilde = Ur'*X2*Vr/Sr;
[W,D] = eig(Atilde);%DMD mode
Phi = X2*Vr/Sr*W;

lambda= diag(D); 
omega = log(lambda)/(dt);%eigenvalue of the solution

%%
x1=X1(:,1);%the sloution at  t=0
b=Phi\x1;%how much time to project
%Time dynamic
time_dynamics=zeros(r,L);
for iter = 1:L
time_dynamics(:,iter) =(b.*exp(omega*t_save(iter)));%make full matrix space and time fot updating t
end
X_dmd = Phi*time_dynamics;


save('SWE_DMD.mat','X_dmd','x','y','H','dt','t_save','plot_height_range')

