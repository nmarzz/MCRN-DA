clear; clc;
% used to create consistent interesting intitial conditions 
% (sets rng seed)
rng(1331);

%% Collect data on full lorenz96 model
[t,y] = ode45(@lorenz96,[0,10],rand(40,1));

X = y';
[U,S,V] = svd(X,'econ');

%% Find and plot singular values and importance
sig=diag(S);
figure(1)
plot(sig,'ko','Linewidth',(1.5)),grid on
xlabel('k')
ylabel('Singular value, \sigma_k')
title('Standard plot of singular values')
 
figure(2)
semilogy(diag(S),'bo','LineWidth',1.5), grid on
xlabel('k')
ylabel('Semilogy of diag(S)')
hold off
title('log plot of singular values')

%% Find how many s values are needed to achieve tol% of the information
Tol=0.5;
cdS =cumsum(sig.^2)./sum(sig.^2);% cumulative 
r =find(cdS>Tol, 1 ); 
figure(3)
plot(cdS,'ko','LineWidth',1.2),grid on
xlabel('k')
ylabel('Cumulative')

%% Truncate matrix
U_r=U(:,1:r);S_r=S(1:r,1:r);V_r=V(:,1:r)';%Truncate U,S,V using the rank r
X_r=U_r*S_r*V_r; %Truncated matrix

%% Try to visualize the POD effect (difficult to do with 40 dimensions)
% for i = 1:length(t) 
%     hold on;
%     plot(t(i),y(i,:),'o');
%     pause(0.05)
% end

















