%% Initialization
clear all;clc;
rng(1331);
% SWE preamble
load('swerun.mat');
F = @(t,x) formod(t,x,dt,pars);
Built_Model= x_save;
N =length(Built_Model);
IC = x_ics;

%% Type of particle filter
% Use of standard PF or OP-PF (iOPPF=0 => standard PF, iOPPF=1 => OP-PF)
iOPPF=0;

%% Projection_type(0 = no projection, 1 POD, 2 DMD, 3 AUS)
PhysicalProjection = 1;
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
Numsteps =Numsteps/2;
% Numsteps = Numsteps/2;
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

gen_ics = x_ics;
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

Q=V'*Q*V; %with projection
Q1=Q;

x=V'*u;
%p_physical = 10;
for i=1:Numsteps
    % Estimate the truth
    estimate(:,i) = x*wt;
    
    % Get projection of the Data Model
    [U] = projectionToggle_data(DataProjection,Ur_data,p_data);
    % Update noise covariance
    UPinvH = U(:,1:inth:end)'; % Multiplying by PinvH is equivalent to removing every inth row when H = H(1:inth:end,:)
    R = UPinvH * Rfixed * UPinvH';
    Rinv = 1/R;
    
    if mod(i,ObsMult)==0
        %At observation times, Update weights via likelihood, add noise
        HV=Hx(V,inth);  % with our assumptions left multipling
        if (iOPPF==0) % Standard Particle Filter
            x = x + mvnrnd(pzeros_physical,Q1,L)';
            Innov=repmat(UPinvH*y(:,i),1,L)-HV*x;
            xx=(1/bet)*Innov;
        else % IOPPF ==1, Optimal proposal PF
            Innov=repmat(UPinvH*y(:,i),1,L)-HV*x;%proj
            Qpinv = (1/alph)*eye(p_physical)+(1/bet)*(HV)'*(HV);
            Qp = pinv(Qpinv);
            x = x + Qp*HV'*Rinv*Innov + mvnrnd(pzeros_physical,Qp,L)';
            yy =(1/bet)*Innov;
            Y = (1/bet)*(HV);
            ww = ((1/alph)+(1/bet)*(HV)'*(HV));
            www=((1/bet)*(HV)'*Innov);
            ww=ww\www;
            xx = (1/bet)*(Innov-(HV)*ww);
            %
            
        end
        
        % Reweight
        Tdiag = diag(Innov'*xx);
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
        %         logw = -(1/2)*1\(alph+bet)+logw;
        wt = exp(logw);
        
        
        %Resampling (with resamp.m that I provided or using the pseudo code in Peter Jan ,... paper)
        [wt,x,NRS] = resamp(wt,x,0.5);
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
    diff_proj= V*(V'* truth(:,i) - estimate(:,i));
    RMSE_orig = sqrt(diff_orig'*diff_orig/N)
    RMSE_proj = sqrt(diff_proj'*diff_proj/N)
    MAE_orig = (sum(abs(diff_orig)))/N;
    RMSEave_orig = RMSEave_orig + RMSE_orig;
    RMSEave_proj = RMSEave_proj + RMSE_proj;
    
    % Save to plot
    if mod(i,ObsMult)==0
        %Save RMSE values
        Time(iRMSE)=t;
        RMSEsave(iRMSE)=RMSE_orig;
        RMSEsave_proj(iRMSE)=RMSE_proj;
        iRMSE = iRMSE+1;
    end
    
    t = t+h;
end
%
figure
plot(Time,RMSEsave, 'r-', 'LineWidth', 1.5)
grid on
hold on;
plot(Time,RMSEsave_proj,'b-','LineWidth', 1.5)
%title('The Root Mean-Squared Error')
xlabel('Time')
ylabel('RMSE')
% xticklabels(xticks/dt)
ylim([0 20])
legend('RMSE Original','RMSE Projected')

RMSEave_orig = RMSEave_orig/Numsteps
RMSEave_proj = RMSEave_proj/Numsteps
ResampPercent = ObsMult*Resamps/Numsteps
