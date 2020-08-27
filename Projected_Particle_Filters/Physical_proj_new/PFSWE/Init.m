function [M,IC,w,R,Rinv,Q,Omega,ICcov,Lones,Mzeros,Nzeros] = Init(IC,N,inth,L,epsR,epsQ,epsOmega,epsIC)

%M = Dimension of observation space.
M = ceil(N/inth);

%Init for Weights
w = zeros(L,1);
w(:)=1/L;

%Init for covariance matrices
%R as the observation error covariance
R = epsR;
Rinv = 1/R;

%Sig as model error covariance
Q = epsQ;

%Omega covariance for Resampling
Omega = epsOmega;

%IC covariance
ICcov = epsIC;

Lones = ones(L,1);
Mzeros=zeros(M,1);

