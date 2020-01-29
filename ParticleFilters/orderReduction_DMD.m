
function [Phi,lambda1] = orderReduction_DMD(r, modeloutput)
X1 = modeloutput(:,1:end-1);%The data matrix
X2 = modeloutput(:,2:end); %shifted data matrix

[U, S, V] = svd(X1,'econ');
r = min(r, size(U,2));
U_r = U(:, 1:r); % truncate to rank -r
S_r = S(1:r, 1:r);
V_r = V(:, 1:r);
Atilde = U_r' * X2 * V_r / S_r; % low -rank dynamics
[W_r , D] = eig(Atilde);
Phi = X2 * V_r / S_r * W_r; % DMD modes
lambda = diag(D); % discrete -time eigenvalues
[lambda1, lambda_index]=sort(real(lambda),'descend');
Phi=Phi(:,lambda_index);%Sort Phi for high value of the real part
% omega = log(lambda)/dt; % continuous -time eigenvalues
end
