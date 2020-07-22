clear all;close all;clc
load('SWE.mat')
% whos
X=modeloutput;%snapshot matrix
L=length(t_save);%number of snapshots


%%
% [U,S,V] = svd(X,'econ');
% figure(3)
% tt= tiledlayout(2, 1);
% nexttile
% plot(real(U(:,1:2)));
% nexttile
% plot(real(V(:,1:2)));
% nexttile([1 2])
% plot((diag(S).^2)/sum(diag(S).^2),'ko','LineWidth',1.2),grid on

%% DMD method
X1 = X(:,1:end-1);%The data matrix
X2 = X(:,2:end); %shifted data matrix

%% Step1:SVD to the data matrix
[U1,S1,V1] = svd(X1, 'econ');
r=60;%number of singular values that we will use

Ur=U1(:,1:r);Sr=S1(1:r,1:r);Vr=V1(:,1:r);%trancated matrices
Atilde = Ur'*X2*Vr/Sr;
[W,D] = eig(Atilde);%DMD mode
Phi = X2*Vr/Sr*W;
lambda= diag(D);
omega = log(lambda)/(dt);%eigenvalue of the solution
%%
x1=X1(:,1);%the sloution at  t=0
b=Phi\x1;%how much time to project
[B,I]=sort(abs(b),'descend');
b=b(I);
Phi=Phi(:,I);
r_small=20;
omega=omega(I);
% b(r_small:end)=0;
tol=1E-8;
if abs(b(r_small) - conj(b(r_small + 1))) < tol
    b(r_small:end)=0;
else
    b(r_small+1:end)=0;
    
end
%Time dynamic
time_dynamics=zeros(r,L);
for iter = 1:L
    time_dynamics(:,iter) =(b.*exp(omega*t_save(iter)));%make full matrix space and time fot updating t
end
X_dmd = Phi*time_dynamics;

save('SWE_DMD_r20.mat','X_dmd','x','y','H','dt','t_save','-v7.3')


