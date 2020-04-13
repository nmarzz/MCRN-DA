function [M,H,PinvH] = new_Init(N,inth)
%Linear Observation operator, every inth variable
Heye=eye(N,N);
H=Heye(1:inth:end,:);
<<<<<<< Updated upstream
% H=Heye(1:inth:end,:)*q;
=======
>>>>>>> Stashed changes
%M = Dimension of observation space.
[M,~]=size(H)
%yvars=linspace(1,M,M);
PinvH=pinv(H);


