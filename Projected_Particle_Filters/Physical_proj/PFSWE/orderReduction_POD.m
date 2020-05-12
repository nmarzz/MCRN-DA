% POD Order Reduction
function [Ur] = orderReduction_POD(r, modeloutput)
[U,S,V] = svd(modeloutput,'econ');
sig=diag(S);
% cdS =cumsum(sig.^2)./sum(sig.^2);% cumulative 
% r = find(cdS>tol, 1 ); 
Ur = U(:,1:r);  
 
% sig=diag(S);
% figure(1)
% plot(sig,'ko','Linewidth',(1.5)),grid on
% xlabel('k')
% ylabel('Singular value, \sigma_k')
% title('Standard plot of singular values')
% May not be necessary
% Sr = S(1:r,1:r);  
% Vr = V(:,1:r)'; % Truncate U,S,V using the rank r
% Xr = Ur*Sr*Vr; % Truncated matrix
end

