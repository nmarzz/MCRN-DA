function [M,IC,w,R,Rinv,Q,Omega,ICcov,Lones,Mzeros,Nzeros] = Init_simp(IC,N,inth,L,epsR,epsQ,epsOmega,epsIC,minidx,maxidx)

%M = Dimension of observation space.
M = size(minidx:inth:maxidx,2);

%Init for Weights
w = zeros(L,1);
w(:)=1/L;

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