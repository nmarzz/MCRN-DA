% POD Order Reduction
function [Ur] = orderReduction_POD(tol, modeloutput)
[U,S] = svd(modeloutput,'econ');
sig=diag(S);
cdS =cumsum(sig.^2)./sum(sig.^2);% cumulative 
r = find(cdS>tol, 1 ); 
Ur = U(:,1:r);  

% May not be necessary
% Sr = S(1:r,1:r);  
% Vr = V(:,1:r)'; % Truncate U,S,V using the rank r
% Xr = Ur*Sr*Vr; % Truncated matrix
end

