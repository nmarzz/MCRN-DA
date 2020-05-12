% POD Order Reduction
function [Ur] = orderReduction_POD(r, modeloutput)
[U] = svd(modeloutput,'econ');
% sig=diag(S);
% cdS =cumsum(sig.^2)./sum(sig.^2);% cumulative 
% r = find(cdS>tol, 1 ); 
Ur = U(:,1:r);  
% Sr = S(1:r,1:r);  
% Vr = V(:,1:r)'; % Truncate U,S,V using the rank r
% Xr = Ur*Sr*Vr; % Truncated matrix
% figure(1)
% plot(sig,'ko','Linewidth',(1.5)),grid on
% xlabel('k')
% ylabel('Singular value, \sigma_k')
% title('Standard plot of singular values')
% 
% figure(2)
% semilogy(diag(S),'bo','LineWidth',1.5), grid on
% xlabel('i')
% ylabel('Semilogy of diag(S)')
% hold off
% title('log plot of singular values')
%  
% figure(3)
% plot(1-cdS,'ko','LineWidth',1.2),grid on
% xlabel('i')
% ylabel('Cumulative')
end
