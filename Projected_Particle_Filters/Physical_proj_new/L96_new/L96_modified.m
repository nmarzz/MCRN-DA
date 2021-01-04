function [Time,RMSEsave, RMSEsave_proj, ResampPercent]=...
    L96_modified(numModes_physical,epsQ, epsR,tolerance_physical,PhysicalProjection,DataProjection,ObsMult,N, F)
% rng(1331);
RHS = @(T,x)FLor95(T,x,F); %Physical model
% N =200; % N:Original model dimension
% Build Model (via ODE45)
dt=1.E-2; % Model output time step
ModelSteps = 50000; % Number of time steps in building model
% ModelSteps = 100000; % Number of time steps in building model
T=ModelSteps*dt;
Built_Model= buildModel(N,RHS,ModelSteps,T);
Built_Model= Built_Model';
IC = zeros(N,1);
IC(1)=1; % Particle ICs
iOPPF=1;

%% Projection_type(0 = no projection, 1 POD, 2 DMD, 3 AUS)
% PhysicalProjection =1;
% DataProjection =1;
% tolerance_physical =5; % POD_modes
tolerance_data =5; % POD_modes
% numModes_physical =6;% DMD_modes, for physical
numModes_data =6; % DMD_modes, for data

model_output = Built_Model;

[Ur_physical,p_physical,pzeros_physical] = ...
    Projection_physical_type(PhysicalProjection ,numModes_physical,tolerance_physical,N,model_output,dt);
[Ur_data,p_data,pzeros_data] = Projection_data_type(DataProjection ,numModes_data,tolerance_data,N,model_output,dt);

%% Particle Filter Information
L=20;%Number of particles
ResampCutoff = 0.5;
% Number of computational steps and step size
% ObsMult=5; % % Observe and every ObsMult steps(10 with F=3.5, 5 with F=8)
h = dt;
NumstepsBig=size(Built_Model,2);
Numsteps=T/h;
alpha=0.99;
% alph=1e-1;
% epsR =0.01;
% epsQ = alph;
epsIC =epsQ;
%%
% IC Variance
epsOmega =1e-2;

inth=1;
[M,IC,wt,R,Rinv,Q,Omega,ICcov,Lones,Mzeros]  = Init_simp_modified(IC,N,inth,L,epsR,epsQ,epsOmega,epsIC);
%Add noise N(0,ICcov) to ICs to form different particles
u = repmat(IC,1,L) + mvnrnd(zeros(1,N),ICcov*ones(1,N),L)'; % Noise for IC
%% Generate observations from "Truth"
y=zeros(size(IC(1:inth:N),1),NumstepsBig);
gen_ics = IC;
t = 0;
for i = 1:NumstepsBig
    truth(:,i) = gen_ics;
    gen_ics =dp4(RHS,t,gen_ics,h);
    if mod(i,ObsMult) == 0
        y(:,i) = truth(1:inth:N,i);
    end
    t = t + h;
end
y = y + mvnrnd(Mzeros',R*ones(1,M),NumstepsBig)'; % Add noise to observations
%% Get projection matrices
[V] =projectionToggle_Physical(PhysicalProjection,N,Ur_physical,p_physical);
[U] =projectionToggle_data(DataProjection,N,Ur_data,p_data);
Rfixed = R;
%%

%Initial time and time step
t=0;
t0=t;
Resamps=0;
RMSEave_orig=0;
RMSEave_proj=0;
% RMSEave_relRMSE=0;
iRMSE=1;
% XC_save_ave=0;
% ESSsave=0;
Q=V'*Q*V; %with projection
Qnew=Q;

x=V'*u;

%Added by EVV
HV=Hx(V,inth,1,N);  % with our assumptions left multipling
UPinvH=Hx(U,inth,1,N)';

Qpfixed = inv(inv(Qnew)+HV'*inv(Rfixed)*HV);
QpHRinv = Qpfixed*HV'*inv(Rfixed);
for i=1:Numsteps
    t;
    if mod(i,ObsMult)==0
        %At observation times, Update particles and weights via likelihood, add noise
        if (iOPPF==0) % Standard Particle Filter
            %Particle Update
            if PhysicalProjection == 0
                [sizepzp,~] = size(pzeros_physical);
                x = x + mvnrnd(pzeros_physical',Qnew*ones(1,sizepzp),L)';
            else
                x = x + mvnrnd(pzeros_physical,Qnew,L)';
            end
            
            %Prepare for weight update
            Innov=repmat(UPinvH*y(:,i),1,L)-UPinvH*HV*x;
            % Update observation covariance
            R = UPinvH * Rfixed * UPinvH';
            RinvtInno = R\Innov;
        else % IOPPF ==1, Optimal proposal PF
            %Particle Update
            Innov1=repmat(y(:,i),1,L)-HV*x;
            if PhysicalProjection == 0
                [sizepzp,~] = size(pzeros_physical);
                x = x + QpHRinv*Innov1 + mvnrnd(pzeros_physical',Qpfixed,L)';
            else
                x = x + QpHRinv*Innov1 + mvnrnd(pzeros_physical,Qpfixed,L)';
            end
            
            %Prepare for weight update
            Innov=UPinvH*Innov1;
            % Update observation covariance
            R = UPinvH * Rfixed * UPinvH';
            Hnq = UPinvH*HV;
            Rnew = R + Hnq*Qnew*Hnq';
            RinvtInno = Rnew\Innov;
        end
        
        % Reweight
        %Tdiag = diag(Innov'*Rinv*Innov);
        Tdiag = diag(Innov'*RinvtInno);
        tempering = 1.2; % Tempering usually a little larger than 1.
        Avg=(max(Tdiag)+min(Tdiag))/2;
        Tdiag = (Tdiag-Avg)/tempering;
        Tdiag = -Tdiag/2;
        
        
        %Take log of the updated weights
        logw = Tdiag + log(wt);
        %To avoid underflow and overflow
        %         logw=min(logw,709);
        %         logw=max(logw,-709);
        %Exponentiate and normalize
        wt = exp(logw);
        wt = wt/sum(wt);
        %Resampling (with resamp.m that I provided or using the pseudo code in Peter Jan ,... paper)
        [wt,x,NRS,ESS] = resamp(wt,x,ResampCutoff);
        Resamps = Resamps + NRS;
        if (NRS==1)
            % Projected Resampling
            x = x + (alpha*(V'*U*U') + (1-alpha)*V')*(mvnrnd(zeros(1,N),epsOmega*ones(1,N),L)');
        end
    end
    
    
    % Estimate the truth
    estimate(:,i) = x*wt;
    x = V'*dp4(RHS,t,V*x,h);
    % Get projection of the Data Model
    % Compare estimate and truth
    diff_orig= truth(:,i) - (V*estimate(:,i));
    diff_proj= V*(V'* truth(:,i) - estimate(:,i));
    Nq=size(V,2);
    RMSE_orig = sqrt(diff_orig'*diff_orig/N);
    RMSE_proj = sqrt(diff_proj'*diff_proj/Nq);
%     xbar=V*estimate(:,i);
%     truth_common=truth(:,i);
%     Ensbar = mean(xbar);
%     Trubar = mean(truth_common);
%     XC_save = (xbar-Ensbar)'*(truth_common-Trubar)/(norm(xbar-Ensbar,2)*norm(truth_common-Trubar));
%     %%
%     Ensbar = mean(V*estimate(:,i));
%     Trubar = mean( truth(:,i));
%     XCproj= (V*estimate(:,i)-Ensbar)'*(truth(:,i)-Trubar)/(norm(V*estimate(:,i)-Ensbar,2)*norm( truth(:,i)-Trubar,2));
    RMSEave_orig = RMSEave_orig + RMSE_orig;
    RMSEave_proj = RMSEave_proj + RMSE_proj;
    
    if mod(i,ObsMult)==0
        %Save RMSE values
        Time(iRMSE)=t;
        RMSEsave(iRMSE)=RMSE_orig;
        RMSEsave_proj(iRMSE)=RMSE_proj;
%         XC_save_ave(iRMSE)=XC_save;
%         XC_save_proj(iRMSE)=XCproj;
%         ESSsave(iRMSE)=ESS;
        iRMSE = iRMSE+1;
    end
    
    t = t+h;
end

RMSEave_orig = RMSEave_orig/Numsteps
RMSEave_proj = RMSEave_proj/Numsteps;
ResampPercent = ObsMult*Resamps/Numsteps*100;
