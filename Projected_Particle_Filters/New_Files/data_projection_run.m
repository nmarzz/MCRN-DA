%% Initialization
clear all;clc;
rng(1331);
Model_Dimension =40;
Fmod = @FLor95;
dt=1.E-2;
Built_Model = buildModel(Model_Dimension,@FLor95,dt);
%% Type of particle filter
% Use of standard PF or OP-PF (iOPPF=0 => standard PF, iOPPF=1 => OP-PF)
iOPPF=0;
proj_data =0;

%% Projection_type(0 = no projection, 1 POD, 2 DMD, 3 AUS)
PhysicalProjection =0;
DataProjection = 1;
tolerance_physical = 0.0001; % POD_modes
tolerance_data = 0.0001; % POD_modes
numModes_physical = 300;% DMD_modes, for physical
numModes_data = 30; % DMD_modes, for data
p_physical = 6;%Rank of projection for physical, in the case PhysicalProjection =0;DataProjection = 0;
p_data = 6; % rank of projection for data
[Ur_physical,p_physical,Nzeros] = ...
    Projection_physical_type(PhysicalProjection ,numModes_physical,tolerance_physical,Model_Dimension,Built_Model,dt,p_physical);
[Ur_data,p_data,Nzeros_data] = Projection_data_type(DataProjection ,numModes_data,tolerance_data,Model_Dimension,Built_Model,dt,p_data);

%% Particle Filter Information
L=2000;%Number of particles
alpha=1;%alpha value for projected resampling
Numsteps = 2000;%Number of time steps
ObsMult=5;%Multiple of the step size for observation time
%ICs for particles
IC = zeros(Model_Dimension,1);
IC(1)=1;

h=1.E-3;%Computational time step
%For diagonal (alpha*I) covariance matrices.
%Observation
epsR = 0.01;
%Model
epsSig = 0.01;
%Resampling
% epsOmega = 2.E-3;
epsOmega =0.0027;
%Initial condition
epsIC = 0.01;
%Observe every inth variable.
inth=1;
%Call Init
[M,H,PinvH,IC,q,LE,w_data,R,Rinv,Sig,Omega,ICcov,Lones,Mzeros] = ...
    Init(Fmod,IC,h,Model_Dimension,inth,Numsteps,p_data,L,epsR,epsSig,epsOmega,epsIC);
Rinvfixed=Rinv;
%Add noise N(0,ICcov) to ICs to form different particles
x_data = zeros(Model_Dimension,L);
Modelzeros = zeros(Model_Dimension,1);
x_data = repmat(IC,1,L) + mvnrnd(Modelzeros,ICcov,L)'; %ICchol*randn(N,L);
y=zeros(M,Numsteps);
t=0;

%% Generate observations from "Truth"
for i = 1:Numsteps
    truth(:,i) = IC;
    if mod(i,ObsMult)==0
        y(:,i)=H*IC + mvnrnd(Mzeros,R,1)'; %Rchol*rand(M,1); % + Noise from N(0,R)
    end
    IC = dp4(Fmod,t,IC,h);
    t = t+h;
end

%%
%Initial time and time step
t=0;
t0=t;
Resamps=0;
RMSEave_orig_data=0;
RMSEave_proj_data=0;
iRMSE=1;
[q_data] = projectionToggle_data(DataProjection,Model_Dimension,Ur_data,p_data); %chooses which q projection we want
%Loop over observation times
[M,H,PinvH] = new_Init(Model_Dimension,inth,q_data);
Sig=q_data'*Sig*q_data;
x_data=q_data'*x_data;
estimate(:,1) = x_data*w_data;
for i=1:Numsteps
    est=estimate(:,i);
    [q_data] = projectionToggle_data(DataProjection,Model_Dimension,Ur_data,p_data); %chooses which q projection we want
    if mod(i,ObsMult)==0
        
        if (iOPPF==0) %Proj-PF
            %H -> Q_n^T P_H where P_H = H^T (H H^T)^{-1} H = H^+ H
            %R -> Q_n^T H^+ R (H^+)^T Q_n where H^+ = H^T (H H^T)^{-1}
            Innov_data = q_data'.*PinvH*(repmat(y(:,i),1,L) - H.*x_data);
            Rinv_data = pinv(q_data'*PinvH*R*PinvH'*q_data);
        else % IOPPF ==1, %Proj-OP-PF
            %H -> Q_n^T P_H where P_H = H^T (H H^T)^{-1} H = H^+ H
            %R -> Q_n^T H^+ R (H^+)^T Q_n where H^+ = H^T (H H^T)^{-1}
            Qpinv_data = inv(Sig) + H'*Rinvfixed*H;
            Qp_data = inv(Qpinv_data);
            Innov_data = repmat(y(:,i),1,L) - H*x_data;
            x_data = x_data + Qp_data*H'*Rinvfixed*Innov_data + mvnrnd(Nzeros_data,Qp_data,L)';
            Innov_data = q_data'*PinvH*Innov_data;
            Rinv_data = pinv(q_data'*PinvH*(R+H*Sig*H')*PinvH'*q_data);
        end
        
        
        Tdiag_data = diag(Innov_data'*Rinv_data*Innov_data);
        tempering = 2; %%%% <<< including new parameter here for visibility. Tempering usually a little larger than 1.
        Avg=(max(Tdiag)+min(Tdiag))/2;
        Tdiag = (Tdiag-Avg)/tempering;
        toEXP_data =(-Tdiag_data/2); %%%% <<<< divided exponent by 2; this is part of the normal distribution
        toEXP_data=max(toEXP_data,-709);%700
        toEXP_data=min(toEXP_data,709);
        LH_data=exp(toEXP_data);
        w_data=LH_data.*w_data;
        %%Normalize weights
        w_data=w_data/(w_data'*Lones);
        %Resampling (with resamp.m that I provided or using the pseudo code in Peter Jan ,... paper)
        [w_data,x_data,NRS_data] = resamp(w_data,x_data,0.5);
        Resamps = Resamps + NRS_data;
        %Update Particles
        
        if (NRS_data==1) %Projected resampling
            x_data = x_data + (alpha*q*q' + (1-alpha))*mvnrnd(Nzeros,Omega,L)'; %Omegachol*randn(N,L);
        end
    end
    %END: At Observation times
    
    
    %Predict, add noise at observation times
    
    % x_data = proj*dp4(Fmod,t,proj*x,h);
    %     x_data = dp4(Fmod,t,x_data,h);
    x_data = q_data'*dp4(Fmod,t,q_data*x_data,h);
    estimate(:,i+1) = x_data*w_data;
    
    diff_orig_data= truth(:,i) - (q_data*estimate(:,i));
    diff_proj= (q_data * q_data'* truth(:,i)) - (q_data*estimate(:,i));
    RMSE_orig_data = sqrt(diff_orig_data'*diff_orig_data/Model_Dimension)
    RMSE_proj_data = sqrt(diff_proj'*diff_proj/Model_Dimension)
    RMSEave_orig_data = RMSEave_orig_data + RMSE_orig_data;
    RMSEave_proj_data = RMSEave_proj_data + RMSE_proj_data;
    
    if mod(i,ObsMult)==0
        %Save RMSE values
        Time(iRMSE)=t;
        RMSEsave(iRMSE)=RMSE_orig_data;
        RMSEsave_proj(iRMSE)=RMSE_proj_data;
        iRMSE = iRMSE+1;
    end
    t = t+h;
end


figure(4)
plot(Time,RMSEsave,'.b');
hold on;
plot(Time,RMSEsave_proj,'.r')
legend('RMSE Original','RMSE Projected')
RMSEave_orig_data = RMSEave_orig_data/Numsteps
RMSEave_proj_data = RMSEave_proj_data/Numsteps
ResampPercent = ObsMult*Resamps/Numsteps


LE = LE/(t-t0)