function [Phi] = orderReduction_DMD(r, modeloutput,dt) % Mar 4, 2020 comments by Marko

X1 = modeloutput(:,1:end-1);%The data matrix
X2 = modeloutput(:,2:end); %shifted data matrix

[U, S, V] = svd(X1,'econ');

figure()
semilogy(diag(S))

% Marko: choice of the following parameters is really not set in stone, compare
Rsmall = min(r, size(U,2)); % order-reduction parameter (modeling choice)
Rlarge = 350 % numerical stabilization (numerical choice)


U_r = U(:, 1:Rlarge); % truncate to rank -r
S_r = S(1:Rlarge, 1:Rlarge);
V_r = V(:, 1:Rlarge);

Atilde = U_r' * X2 * V_r / S_r; % low -rank dynamics
[W_r , D] = eig(Atilde);
Phi = X2 * V_r / S_r * W_r; % DMD modes
lambda = diag(D); % discrete -time eigenvalues

omega = log(lambda)/dt; % continuous -time eigenvalues

% THE TRUE SORTING HAPPENS HERE
[~, omega_index]=sort(real(omega),'descend');

omega = omega(omega_index); % sort continuous-time eigenvalues
Phi=Phi(:,omega_index);% sort DMD modes

% this is one option: keep only real parts of modes
Phi=real(Phi);

% the other option is to manually go through columns/omega eigenvalues
% and replace pair of conjugate modes with real() , imag() of Phi
% Notice that if you take real and imaginary parts from
% v = a + ib and
% conj(v) = a - ib
% you get
% * (a, b) [if you take from v] or
% * (a, -b) [if you take from conj(v)]
% In both cases, the same space is spanned.

% make sure I don't use same mode twice
CC=round(Phi,6);% rounding
[~,index,~]=unique(CC.','rows');
index_sort=sort(index);

Phi=Phi(:,index_sort);


Phi = Phi(:, 1:Rsmall);
size(Phi)
rank(Phi)

end