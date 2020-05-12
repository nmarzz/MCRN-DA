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


%% Projection_type(0 = no projection, 1 POD, 2 DMD, 3 AUS)
PhysicalProjection =2;
DataProjection = 2;
tolerance_physical = 0.0001; % POD_modes
tolerance_data = 0.0001; % POD_modes
numModes_physical = 30;% DMD_modes, for physical
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
inth=2;
%Call Init
[M,H,PinvH,IC,q,LE,w,R,Rinv,Sig,Omega,ICcov,Lones,Mzeros] = ...
    Init(Fmod,IC,h,Model_Dimension,inth,Numsteps,p_physical,L,epsR,epsSig,epsOmega,epsIC);
% Rinvfixed=Rinv;
Rfixed=R;
%Add noise N(0,ICcov) to ICs to form different particles
x = zeros(Model_Dimension,L);
Modelzeros = zeros(Model_Dimension,1);
x = repmat(IC,1,L) + mvnrnd(Modelzeros,ICcov,L)'; %ICchol*randn(N,L);
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
RMSEave_orig=0;
RMSEave_proj=0;
iRMSE=1;
[q_physical] =projectionToggle_Physical(PhysicalProjection,Model_Dimension,Ur_physical,p_physical);
%Loop over observation times
[M,H,PinvH] = new_Init(Model_Dimension,inth);
Sig=q_physical'*Sig*q_physical;
x=q_physical'*x;
estimate(:,1) = x*w;

for i=1:Numsteps
    % Perform projection of the Data Model
    est=estimate(:,i);
    [q_data] = projectionToggle_data(DataProjection,Model_Dimension,Ur_data,p_data); %chooses which q projection we want
    R=q_data'*PinvH*Rfixed*PinvH'*q_data;%udataing R
    if mod(i,ObsMult)==0
        %At observation times, Update weights via likelihood
        if (iOPPF==0)%Standard Particle Filter(Physical projection)
            %Add noise only at observation times
            x = x + mvnrnd(Nzeros,Sig,L)';
            %             Innov = repmat(y(:,i),1,L) - H*x;
            % Innov = q_data'*PinvH*(repmat(y(:,i),1,L) - H*x);%Proj_data_model
            Innov = q_data'*PinvH*(repmat(y(:,i),1,L) - H*q_physical*x);%Proj_data_model and Physical
        else % IOPPF ==1, Optimal proposal PF
            Qpinv = inv(Sig) + H'*Rinvfixed*H;
            Qp = inv(Qpinv);
            Innov = repmat(y(:,i),1,L) - H*x;
            %             Innov = q_data'*PinvH*(repmat(y(:,i),1,L) - H*q_physical*x);
            x = x + Qp*H'*Rinv*Innov + mvnrnd(Nzeros,Qp,L)';
            Rinv = inv(R + H*Sig*H');
%             Innov = repmat(y(:,i),1,L) - H*x;
            Hnq = q_data'*Hcross*H*q_physical; % H has already been multiplied by q_physical in new_init            
            Innov=q_data'*Hcross*repmat(y(:,i),1,L)-Hnq*x;%try to code what Erik wrote on slack
        else % IOPPF ==1, Optimal proposal PF
            Hnq = q_data'*Hcross*H*q_physical; % H has already been multiplied by q_physical in new_init  
            Qpinv = inv(Sig) + Hnq'*Rinv*Hnq;
            Qp = inv(Qpinv);
            Innov=q_data'*Hcross*repmat(y(:,i),1,L)-Hnq*x;
            x = x + Qp*Hnq'*Rinv*Innov + mvnrnd(Nzeros,Qp,L)';
            Rinv = inv(R + Hnq*Sig*Hnq');

        end
        Tdiag = diag(Innov'*Rinv*Innov);
        tempering = 2; % including new parameter here for visibility. Tempering usually a little larger than 1.
        Avg=(max(Tdiag)+min(Tdiag))/2;
        Tdiag = (Tdiag-Avg)/tempering;
        %         % NEW CODE: Re weight while avoiding taking large exponentials -
        %         % avoids NAN more often
        Tdiag = -Tdiag/2;
        logw = Tdiag + log(w);
        
        % Identity used to redo normalization: log(sum(a)) = log(a_0) + log(1 + sum exp(log(w0(2:end)) - log(w0(1))))
        % Re-order weights to ensure we take the smallest possible exponents
        [~,idx] = min(abs(logw-((max(logw) - min(logw))/2))); % find index of weight closest to middle value
        logw([1 idx]) = logw([idx 1]);
        x(:,[1 idx]) = x(:,[idx 1]);
        
        toEXP = logw(2:end) - logw(1);
        toEXP=min(toEXP,709);
        toEXP=max(toEXP,-709);
        toSum = exp(toEXP);
        normalizer = logw(1) + log1p(sum(toSum));
        logw = logw - normalizer;
        w = exp(logw);
        
        % End new code
        
        
        %%%LEGACY CODE - normalizing and reweighting
        %         toEXP =(-Tdiag/2); %%%% <<<< divided exponent by 2; this is part of the normal distribution
        %         toEXP=max(toEXP,-709);%700
        %         toEXP=min(toEXP,709);
        %         LH=exp(toEXP);
        %         w=LH.*w;
        %         %%Normalize weights
        %         w=w/(w'*Lones);
        %Resampling (with resamp.m that I provided or using the pseudo code in Peter Jan ,... paper)
        [w,x,NRS] = resamp(w,x,0.5);
        Resamps = Resamps + NRS;
        %Update Particles
        %Note: This can be modified to implement the projected resampling as part of implementation of PROJ-PF:
        %Replace Sigchol*randn(N,L) with (alpha*Q_n*Q_n^T + (1-alpha)I)*Sigchol*randn(N,L)
        if (NRS==1)
            %Standard resampling
            x = x + mvnrnd(Nzeros,Sig,L)'; %Sigchol*randn(N,L);
        end
        %END: At Observation times
    end
    
    %Predict, add noise at observation times
    %     x = proj*dp4(Fmod,t,proj*x,h);
    x = q_physical'*dp4(Fmod,t,q_physical*x,h);
    
    estimate(:,i+1) = x*w;
    
    diff_orig= truth(:,i) - (q_physical*estimate(:,i));
    diff_proj= (q_physical * q_physical'* truth(:,i)) - (q_physical*estimate(:,i));
    RMSE_orig = sqrt(diff_orig'*diff_orig/Model_Dimension)
    RMSE_proj = sqrt(diff_proj'*diff_proj/Model_Dimension)
    RMSEave_orig = RMSEave_orig + RMSE_orig;
    RMSEave_proj = RMSEave_proj + RMSE_proj;
    
    if mod(i,ObsMult)==0
        %Save RMSE values
        Time(iRMSE)=t;
        RMSEsave(iRMSE)=RMSE_orig;
        RMSEsave_proj(iRMSE)=RMSE_proj;
        iRMSE = iRMSE+1;
    end
    t = t+h;
end
figure(4)
plot(Time,RMSEsave,'.b');
hold on;
plot(Time,RMSEsave_proj,'.r')
legend('RMSE Original','RMSE Projected')
RMSEave_orig = RMSEave_orig/Numsteps
RMSEave_proj = RMSEave_proj/Numsteps
ResampPercent = ObsMult*Resamps/Numsteps


