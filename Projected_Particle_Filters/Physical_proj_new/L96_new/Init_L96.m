function [M,H,PinvH,IC,q,LE,w,R,Rinv,Q,Omega,ICcov,Lones,Mzeros,Nzeros] = Init_L96(F,IC,h,N,inth,Numsteps,p,L,epsR,epsQ,epsOmega,epsIC)

%Linear Observation operator, every inth variable
Heye=eye(N,N);
H=Heye(1:inth:end,:);
%M = Dimension of observation space.
[M,~]=size(H)
%yvars=linspace(1,M,M);
PinvH=pinv(H);

%Init for LEs and projections
q=orth(randn(N,p));
LE=zeros(p,1);

%Spin up
t=0;
%for i = 1:Numsteps*2
for i = 1:5000
t = t+h;
[q,LE] = getausproj(N,p,F,t,IC,h,q,LE);
IC = dp4(F,t,IC,h);
end

LE=zeros(p,1);

%Init for Weights
w = zeros(L,1);
w(:)=1/L;

%Init for covariance matrices
%R as the observation error covariance
R = epsR*eye(M);
Rinv = inv(R);

%Sig as model error covariance
Q = epsQ*eye(N);

%Omega covariance for Resampling
Omega = epsOmega*eye(N);

%IC covariance
ICcov = epsIC*eye(N);

Lones = ones(L,1);
Mzeros=zeros(M,1);

