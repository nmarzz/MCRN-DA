clearvars; close all;clc;
load('SWE_run_4days.mat');
parsanim.H = zeros(pars.nx, pars.ny);
parsanim.x=(0:pars.nx-1).*pars.dx; % Zonal distance coordinate (m)
parsanim.y=(0:pars.ny-1).*pars.dy; % Meridional distance coordinate (m)
F = @(t,x) formod(t,x,dt,pars);
Built_Model= x_save;
N =length(Built_Model);
IC = Built_Model(:,(end-1)/2);
%% Type of particle filter
% Use of standard PF or OP-PF (iOPPF=0 => standard PF, iOPPF=1 => OP-PF)
iOPPF=1;
scenario = 2;
[minidx,maxidx] = getScenarioIndex(scenario,N);
%% Projection_type(0 = no projection, 1 POD, 2 DMD, 3 AUS)
PhysicalProjection =0;
DataProjection =0;
numModes_physical =11;% DMD_modes, for physical
tolerance_physical =20; % POD_modes
tolerance_data =10; % POD_modes
numModes_data =11; % DMD_modes, for data
model_output = Built_Model(:,(end-1)/2:(end-1)*3/4);
[Ur_physical,p_physical,pzeros_physical] = ...
    Projection_physical_type(PhysicalProjection ,numModes_physical,tolerance_physical,N,model_output,dt);
[Ur_data,p_data,pzeros_data] = Projection_data_type(DataProjection ,numModes_data,tolerance_data,N,model_output,dt);

%% Particle Filter Information
ResampCutoff = 0.5;
% Number of computational steps and step size
ObsMult=60; % Observe and every ObsMult steps
h = dt;
Numsteps=1441;
NumstepsBig=size(Built_Model,2);
epsOmega =0.0000001; %For inth = 1
%Observe every inth variable.
inth=100;%observations run for 10,100,1000
epsQ=1E-1; %Run for epsQ=1
epsR=1E-2;
L=5; % number of particles



%% For PF-OP
alpha =.99;%alpha value for projected resampling
epsIC = epsQ;
%%
%Call Init
[M,IC,wt,R,Rinv,Q,Omega,ICcov,Lones,Mzeros] = Init_simp(IC,N,inth,L,epsR,epsQ,epsOmega,epsIC,minidx,maxidx);
%Add noise N(0,ICcov) to ICs to form different particles
u = repmat(IC,1,L) + mvnrnd(zeros(1,N),ICcov*ones(1,N),L)'; % Noise for IC
%% Generate observations from "Truth"
y=zeros(size(IC(minidx:inth:maxidx),1),NumstepsBig);
gen_ics = IC;
t = 0;
for i = 1:NumstepsBig
    truth(:,i) = gen_ics;
    gen_ics = formod(t,gen_ics,dt,pars);
    if mod(i,ObsMult) == 0
        y(:,i) = truth(minidx:inth:maxidx,i);
    end
    t = t + h;
end
y = y + mvnrnd(Mzeros',R*ones(1,M),NumstepsBig)'; % Add noise to observations
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
% RMSEave_proj=0;
% RMSEave_relRMSE=0;
iRMSE=1;
% % XC_save_ave=0;
ESSsave=0;
Q=V'*Q*V; %with projection
Qnew=Q;

x=V'*u;

%Added by EVV
HV=Hx(V,inth,minidx,maxidx);  % with our assumptions left multipling
UPinvH=Hx(U,inth,minidx,maxidx)';

Qpfixed = inv(inv(Qnew)+HV'*inv(Rfixed)*HV);
QpHRinv = Qpfixed*HV'*inv(Rfixed);
% diff_plot=[];
% diff_proj_plot=[];
% XC_save=[];
% e.g. lowerbound( [0.1, 0.2, 0.3], 0.15 ) = [0.15, 0.2, 0.3]
% lowerbound = @(x,lb) max( ...
%     cat(length(size(x))+1, ... % concatenate along the "last" dimension (e.g.,
%     x, ... % input vector
%     repmat(lb, size(x)) ... % vector of lower bounds in the shape of input x
%     ), ...
%     [], length(size(x))+1 ) ;% and then take the max of both.
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
                % x = x + QpHRinv*Innov1 + mvnrnd(pzeros_physical',(Qpfixed*ones(1,sizepzp)')',L)';
                x = x + QpHRinv*Innov1 + mvnrnd(pzeros_physical',(Qpfixed*ones(1,sizepzp)')',L)';
            else
                x = x + QpHRinv*Innov1 + mvnrnd(pzeros_physical,Qpfixed,L)';
            end
            
            %Prepare for weight update
            Innov=UPinvH*Innov1;
            % Update observation covariance
            R = UPinvH * Rfixed * UPinvH';
            Hnq = UPinvH*HV;
            Rnew = R + Hnq*Qnew*Hnq';
%             RinvtInno = Rnew\Innov;
            RinvtInno = pinv(full(Rnew))*Innov;
%             any(isnan(RinvtInno))
%             RinvtInno = RinvtInno*Innov;
            
        end
        
        % Reweight
        %         Tdiag = diag(Innov'*Rinv*Innov);
        Tdiag = diag(Innov'*RinvtInno);
        any(isnan(Tdiag));
        tempering = 1.2; % Tempering usually a little larger than 1.
        Avg=(max(Tdiag)+min(Tdiag))/2;
        Tdiag = (Tdiag-Avg)/tempering;
        Tdiag = -Tdiag/2;
        
        
        %Take log of the updated weights
        logw = Tdiag + log(wt);
        
        %% normalize the weights
        maxlogw = max(logw);
        
        logw = normlogw( logw ); % normalization w/o overflow
        
        wt = exp(logw);
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
    
    % Get projection of the Data Model
    for k = 1:size(x,2) % formod isn't vectorized
        % propogate particles
        x(:,k) = V'*formod(t,V*x(:,k),dt,pars);
    end
    % Compare estimate and truth
    diff_orig= truth(:,i) - (V*estimate(:,i));
    %     diff_proj= V*(V'* truth(:,i) - estimate(:,i));
    RMSE_orig = sqrt(diff_orig'*diff_orig/N);
    
    
    %     Nq=size(V,2);
    %     RMSE_proj = sqrt(diff_proj'*diff_proj/Nq);
    %xbar=V*estimate(:,i);
    %truth_common=truth(:,i);
    %Ensbar = mean(xbar);
    %Trubar = mean(truth_common);
    %XC_save = (xbar-Ensbar)'*(truth_common-Trubar)/(norm(xbar-Ensbar,2)*norm(truth_common-Trubar));
    %%
    %Ensbar = mean(V*estimate(:,i));
    %Trubar = mean( truth(:,i));
    %XCproj= (V*estimate(:,i)-Ensbar)'*(truth(:,i)-Trubar)/(norm(V*estimate(:,i)-Ensbar,2)*norm( truth(:,i)-Trubar,2));
    
    %     xbar=V*estimate(:,i);
    %     truth_common=V*V'*truth(:,i);
    %     Ensbar = mean(xbar);
    %     Trubar = mean(truth_common);
    %     XCproj = (xbar-Ensbar)'*(truth_common-Trubar)/(norm(xbar-Ensbar,2)*norm(truth_common-Trubar));
    RMSEave_orig = RMSEave_orig + RMSE_orig;
    
    disp("RMSEave_orig = " + RMSEave_orig + ...
        "  weight_min = " + min(wt) + ...
        "  weight_max = " + max(wt) );
    
    
    if mod(i,ObsMult)==0
        %Save RMSE values
        Time(iRMSE)=t;
        RMSEsave(iRMSE)=RMSE_orig;
        %         RMSEsave_proj(iRMSE)=RMSE_proj;
        %         XC_save_ave(iRMSE)=XC_save;
        %         XC_save_proj(iRMSE)=XCproj;
        ESSsave(iRMSE)=ESS;
        iRMSE = iRMSE+1;
    end
    
    t = t+h;
end

RMSEave_orig = RMSEave_orig/Numsteps
% RMSEave_proj = RMSEave_proj/Numsteps;
ResampPercent = ObsMult*Resamps/Numsteps*100;
