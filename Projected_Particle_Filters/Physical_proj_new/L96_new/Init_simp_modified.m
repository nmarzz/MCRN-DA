function [M,IC,w,R,Rinv,Q,Omega,ICcov,Lones,Mzeros] = Init_simp_modified(IC,N,inth,L,epsR,epsQ,epsOmega,epsIC,F,h)

%M = Dimension of observation space.
M = size(1:inth:N,2);

%Init for Weights
w = zeros(L,1);
w(:)=1/L;
t=0;
for i=1:1000
IC =dp4(F,t,IC,h);
t=t+h;
end
%Init for covariance matrices
%R as the observation error covariance
R = epsR;
%*eye(N);
%epsR;
Rinv = 1/R;

%Sig as model error covariance
Q = epsQ;
%* eye(N);
%epsQ;

%Omega covariance for Resampling
Omega = epsOmega;

%IC covariance
ICcov = epsIC;

Lones = ones(L,1);
Mzeros=zeros(M,1);