function [Phi]=buildDMD(numModes,model_output,dt)
% DMD
% model_output = model_output';
% [Phi] = orderReduction_DMD(numModes, model_output,dt);% Get DMD
X1 = model_output(:,1:end-1);%The data matrix
X2 = model_output(:,2:end); %shifted data matrix
[U, S, V] = svd(X1,'econ');
%%
%sig=diag(S);
% figure(1)
% plot(sig,'ko','Linewidth',(1.5)),grid on
% xlabel('k')
% ylabel('Singular value, \sigma_k')
% title('Standard plot of singular values')

% figure(2)
% semilogy(diag(S),'bo','LineWidth',1.5), grid on
% xlabel('k')
% ylabel('Semilogy of diag(S)')
% hold off
% title('log plot of singular values')
% cdS =cumsum(sig.^2)./sum(sig.^2);% cumulative
% figure(3)
% plot(cdS,'ko','LineWidth',1.2),grid on
% xlabel('k')
% ylabel('Cumulative')
%%
% Marko: choice of the following parameters is really not set in stone, compare
Rsmall = min(numModes, size(U,2)); % order-reduction parameter (modeling choice)
Rlarge = size(U, 2) - 10; % numerical stabilization (numerical choice)
% Rlarge = rank(modeloutput) - 10;


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

% omega = omega(omega_index); % sort continuous-time eigenvalues
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

Rsmall=min(Rsmall,size(index_sort));
Phi = Phi(:, 1:Rsmall);
size(Phi)
rank(Phi)
end
