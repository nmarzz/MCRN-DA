%% Initialization
clear all;clc;
rng(1331);
% SWE preamble
% load('swerun.mat');
% F = @(t,x) formod(t,x,dt,pars);
% Built_Model= x_save;
% N =length(Built_Model);
% IC = x_ics;
load('SWE_RUN_2days.mat');
F = @(t,x) formod(t,x,dt,pars);
Built_Model= x_save;
N =length(Built_Model);
IC =Built_Model(:,1);
% model=Built_Model(1:2000,1:2000);
% figure(9)
% contourf(model,'LineStyle','none')
%% Type of particle filter
% Use of standard PF or OP-PF (iOPPF=0 => standard PF, iOPPF=1 => OP-PF)
iOPPF=0;

%% Projection_type(0 = no projection, 1 POD, 2 DMD, 3 AUS)
PhysicalProjection = 0;
DataProjection = 0;
tolerance_physical = 9; % POD_modes
tolerance_data = 10; % POD_modes
numModes_physical = 30;% DMD_modes, for physical
numModes_data = 30; % DMD_modes, for data

model_output = Built_Model;
[Ur_physical,p_physical,pzeros_physical] = ...
    Projection_physical_type(PhysicalProjection ,numModes_physical,tolerance_physical,N,model_output,dt);
[Ur_data,p_data,pzeros_data] = Projection_data_type(DataProjection ,numModes_data,tolerance_data,N,model_output,dt);

%% Particle Filter Information
L=200;%Number of particles
alph = 0.01;
bet = 0.01;
alpha = 0;%alpha value for projected resampling
% Number of computational steps and step size
ObsMult=1; % Observe and every ObsMult steps
h = dt/ObsMult;
Numsteps = size(Built_Model,2)*ObsMult;
Numsteps = Numsteps/2;
%Observation Variance
epsR = 0.01;
%Model Variance
epsQ = 0.01;
% IC Variance
epsOmega =0.0027;
%Initial condition
epsIC = 0.01;
%Observe every inth variable.
inth=2;
%Call Init
[M,IC,wt,R,Rinv,Q,Omega,ICcov,Lones,Mzeros] = Init_simp(IC,N,inth,L,epsR,epsQ,epsOmega,epsIC,alph, bet);
%Init(IC,N,inth,L,epsR,epsQ,epsOmega,epsIC);

%Add noise N(0,ICcov) to ICs to form different particles

u = repmat(IC,1,L) + normrnd(0,ICcov,N,L); % Noise for IC

%% Generate observations from "Truth"

gen_ics = IC;
t = 0;
for i = 1:Numsteps
    truth(:,i) = gen_ics;
    gen_ics = formod(t,gen_ics,dt,pars);
    if mod(i,ObsMult) == 0
        y = truth(1:inth:end,i);
    end
    t = t + h;
end
y = y + normrnd(0,R,M,size(Built_Model,2));
%% Get projection matrices
[V] =projectionToggle_Physical(PhysicalProjection,Ur_physical,p_physical);
[U] =projectionToggle_data(DataProjection,Ur_data,p_data);
Rfixed = R;
%%
%Initial time and time step
t=0;
t0=t;
Resamps=0;
RMSEave_orig=0;
RMSEave_proj=0;
iRMSE=1;
Q1=Q*ones(1,N); %no Projection
x=V'*u;
%p_physical = 10;
for i=1:Numsteps
    % Estimate the truth
    estimate(:,i) = x*wt;
    
    % Get projection of the Data Model
    [U] = projectionToggle_data(DataProjection,Ur_data,p_data);
    % Update noise covariance
    UPinvH = U(:,1:inth:end)'; % Multiplying by PinvH is equivalent to removing every inth row when H = H(1:inth:end,:)
    % UPinvH = U' * PinvH = U' * pinv(H);
    R = UPinvH * Rfixed * UPinvH';
    Rinv = 1/R;
    
    if mod(i,ObsMult)==0
        %At observation times, Update weights via likelihood, add noise
        if (iOPPF==0) % Standard Particle Filter
            x = x + mvnrnd(pzeros_physical,Q1,L)';
            Hnq=Hx(V,inth); % with our assumptions left multipling
            Innov=repmat(UPinvH*y(:,i),1,L)-x(1:inth:end,:);
        else % IOPPF ==1, Optimal proposal PF
            Hnq=Hx(V,inth);
            Innov=repmat(y(:,i),1,L)-Hnq*x(1:inth:end,:);
            Qpinv = (1/alph)+(1/bet)*(Hnq*V)'*(Hnq*V);
            Qp = pinv(Qpinv);
            Qp1=Qp*eye(1,N);
            x = x + (alph*bet)/(alph+bet)* Rinv* Hty(Innov,inth,N,L)+ mvnrnd(pzeros_physical,Qp1,L)';
            Rinv = inv(R + Hnq*Q*Hnq');
        end
        
        % Reweight
        Tdiag = diag(Innov'*Rinv*Innov);
        tempering = 2; % Tempering usually a little larger than 1.
        Avg=(max(Tdiag)+min(Tdiag))/2;
        Tdiag = (Tdiag-Avg)/tempering;
        
        Tdiag = -Tdiag/2;
        logw = Tdiag + log(wt);
        
        [~,idx] = min(abs(logw-((max(logw) - min(logw))/2))); % find index of weight closest to middle value
        logw([1 idx]) = logw([idx 1]);
        x(:,[1 idx]) = x(:,[idx 1]);
        toEXP = logw(2:end) - logw(1);
        toEXP=min(toEXP,709);
        toEXP=max(toEXP,-709);
        toSum = exp(toEXP);
        normalizer = logw(1) + log1p(sum(toSum));
        logw = logw - normalizer;
        wt = exp(logw);
        
        
        %Resampling (with resamp.m that I provided or using the pseudo code in Peter Jan ,... paper)
        [wt,x,NRS] = resamp(wt,x,0.3);
        Resamps = Resamps + NRS;
        
        if (NRS==1)
            % Projected Resampling
            x = x + V'*(alpha*(U*U') + (1-alpha)*eye(N,1))*(mvnrnd(0,Omega,L)');% with proj
            
        end
        
    end
    
    for k = 1:size(x,2) % formod isn't vectorized
        % propogate particles
        x(:,k) = V'*formod(t,V*x(:,k),dt,pars);
    end
    % Compare estimate and truth
    diff_orig= truth(:,i) - (V*estimate(:,i));
    RMSE_orig = sqrt(diff_orig'*diff_orig/N);
    %     MAE_orig = (sum(abs(diff_orig)))/N
    RMSEave_orig = RMSEave_orig + RMSE_orig;
    
    % Save to plot
    if mod(i,ObsMult)==0
        %Save RMSE values
        Time(iRMSE)=t;
        RMSEsave(iRMSE)=RMSE_orig;
        iRMSE = iRMSE+1;
    end
    
    t = t+h;
end
%
% epsRR=epsR*ones(1,length(RMSEsave));
% % loyolagreen = 1/255*[0,104,87];
% figure
% plot(Time,RMSEsave, 'b--', 'LineWidth', 1.5)
% grid on
% hold on;
% % plot(Time,RMSEsave_proj,'r--','LineWidth', 1.5)
% % plot(Time,epsRR,':','Color', loyolagreen,'LineWidth', 1.5)
% plot(Time,epsRR,'g-.','LineWidth', 1.5)
% xlabel('Time')
% ylabel('RMSE')
% ylim([0 25])
% title('Identity Projection')
% % title('DMD Projection')
% % title('POD and DMD Projection')
% % legend('RMSE','Observation error','Location', 'Best')
% legend('RMSE','Observation error','Location', 'Best')
% 
% 
RMSEave_orig = RMSEave_orig/Numsteps
ResampPercent = ObsMult*Resamps/Numsteps
