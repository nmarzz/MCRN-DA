function [M,H,PinvH] = new_Init(N,inth,q)
%Linear Observation operator, every inth variable
Heye=eye(N,N);
H=Heye(1:inth:end,:)*q;
%M = Dimension of observation space.
[M,~]=size(H)
%yvars=linspace(1,M,M);
PinvH=pinv(H);


