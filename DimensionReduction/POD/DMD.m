clc;clear all;close all;
load('rotor_oscillator_vf.mat')
whos
%The data have:
%Ux-Horizontal velocity components
%Uy-Vertical velocity components
%x-dim,y-dim for each dt-time step and t is the time
%% Step1,build a matrix such that time x-dir and space y-dir 
L=length(t);%number of snapshots 
m=length(x);%y-dimension 
n=length(y);%x-dimension
k=m*n;%each snapshot has x*y-dimensions
Ux1=zeros(k,L);%store new Ux matrix 
Uy1=zeros(k,L);%store new Uy matrix

for i=1:L
    Ux1(:,i) =reshape(Ux(:,:,i),[k,1]);
    Uy1(:,i) =reshape(Uy(:,:,i),[k,1]);
end
 X=[Ux1;Uy1]; %the time-space matrix 
 
%%
[U,S,V] = svd(X,'econ');
figure(3)
tt= tiledlayout(2, 1);
nexttile
plot(real(U(:,1:2)));
nexttile
plot(real(V(:,1:2)));
% nexttile([1 2])
% plot((diag(S).^2)/sum(diag(S).^2),'ko','LineWidth',1.2),grid on

%% DMD method
X1 = X(:,1:end-1);%The data matrix
X2 = X(:,2:end); %shifted data matrix

%% Step1:SVD to the data matrix
[U1,S1,V1] = svd(X1, 'econ');
r=2;%number of singular values that we will use
Ur=U1(:,1:r);Sr=S1(1:r,1:r);Vr=V1(:,1:r);%trancated matrices
Atilde = Ur'*X2*Vr/Sr;
[W,D] = eig(Atilde);%DMD mode
Phi = X2*Vr/Sr*W;

lambda= diag(D); 
omega = log(lambda)/(dt);%eigenvalue of the solution
%%
figure(4)
tlt = tiledlayout(2, 2);
nexttile
plot(real(U(:,1:r)));
nexttile
plot(real(Phi));
nexttile([1 2])
plot(real(U(:,1:r)),'-.','Linewidth',1.5);
hold on
plot(real(Phi));
%%
x1=X1(:,1);%the sloution at  t=0
b=Phi\x1;%how much time to project
%Time dynamic
time_dynamics=zeros(r,length(t));
for iter = 1:length(t)
time_dynamics(:,iter) =(b.*exp(omega*t(iter)));%make full matrix space and time fot updating t
end
X_dmd = Phi*time_dynamics;
%% Step3: Regenarate the horizontal and vertical velocity components

Ux2=X_dmd(1:k,:);
Uy2=X_dmd(k+1:end,:);

Ux3=reshape(Ux2,n,m,L);
Uy3=reshape(Uy2,n,m,L);

save('velocity_dmd','Ux3','Uy3','x','y','t','dt')

