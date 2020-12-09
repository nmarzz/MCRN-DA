function [Time,RMSEsave, RMSEsave_proj, XCsave, XCprojsave, ESSsave, ResampPercent]=...
    PFL96run_1(tolerance_physical,numModes_physical,phys_proj,data_proj,epsQ,epsR,iOPPF,N)
rng(1331);
F = @FLor95; %Physical model
% N =100; % N:Original model dimension
% Build Model (via ODE45)
dt=1.E-2; % Model output time step
ModelSteps = 50000; % Number of time steps in building model
% ModelSteps = 100000; % Number of time steps in building model
T=ModelSteps*dt;
Built_Model= buildModel(N,F,ModelSteps,T);
model_output = Built_Model';

%% Projection_type(0 = no projection, 1 POD, 2 DMD, 3 AUS)
PhysicalProjection =phys_proj;
DataProjection = data_proj;
tolerance_data = 5; % POD_modes
%numModes_physical = 20;% DMD_modes/AUS_modes, for physical
%numModes_data = 20; % DMD_modes/AUS_modes, for data
numModes_data = 5; % DMD_modes/AUS_modes, for data

[Ur_physical,p_physical,pzeros_physical] = ...
    Projection_physical_type(PhysicalProjection, numModes_physical,tolerance_physical,N,model_output,dt);
[Ur_data,p_data,pzeros_data] = Projection_data_type(DataProjection ,numModes_data,tolerance_data,N,model_output,dt);

%% Particle Filter Information
L=20;%Number of particles
IC = zeros(N,1);
IC(1)=1; % Particle ICs

%alpha=0.35;%alpha value for projected resampling
alpha=0.99;%alpha value for projected resampling
% ResampCutoff=0.3;
ResampCutoff=0.5;

% Number of computational steps and step size
h = dt;
Numsteps=T/h;

ObsMult=5; % Observe and every ObsMult steps(10 with F=3.5, 5 with F=8)

%Observation Variance
epsOmega =1e-2;
epsIC =epsQ;
%inth=2;
inth=1;
%Call Init
if PhysicalProjection == 3
    NumLEs=p_physical;
else
    NumLEs=p_data;
end
[M,H,PinvH,IC,q,LE,w,R,Rinv,Q,Omega,ICcov,Lones,Mzeros] = ...
    Init_L96(F,IC,h,N,inth,ModelSteps,NumLEs,L,epsR,epsQ,epsOmega,epsIC);

%Add noise N(0,ICcov) to ICs to form different particles
Nzeros = zeros(N,1);
u = repmat(IC,1,L) + mvnrnd(Nzeros,ICcov,L)'; % Noise for IC

%% Generate observations from "Truth"
y=zeros(M,Numsteps);
t=0;
for i = 1:Numsteps
    truth(:,i) = IC;
    if mod(i,ObsMult)==0
        y(:,i)=H*IC;
    end
    IC = dp4(F,t,IC,h);
    t = t+h;
end
y = y + mvnrnd(Mzeros,R,Numsteps)'; % Add noise to observations

%% Get projection matrices
if (PhysicalProjection == 3)
    V=q;
else
    [V] =projectionToggle_Physical(PhysicalProjection,N,Ur_physical,p_physical);
    V0=V;
end
if (DataProjection == 3)
    U=q;
else
    [U] =projectionToggle_data(DataProjection,N,Ur_data,p_data);
end
%%

Rfixed = R;

%Initial time and time step
t=0;
t0=t;
Resamps=0;
RMSEave_orig=0;
RMSEave_proj=0;
XC=0;
ESS=0;
XCproj=0;
iRMSE=1;
[M,H,PinvH] = new_Init_L96(N,inth,V);
x=V'*u;

Qfixed=Q;
Q=V'*Q*V;

%Added by EVV
Qnew=V'*Qfixed*V;
Qpfixed = inv(inv(Qnew)+V'*H'*inv(Rfixed)*H*V);
QpHRinv = Qpfixed*V'*H'*inv(Rfixed);


for i=1:Numsteps
    
    if mod(i,ObsMult)==0
        %At observation times, Update weights via likelihood, add noise
        if (iOPPF==0) % Standard Particle Filter
            x = x + mvnrnd(pzeros_physical,Q,L)';
            Hnq = U'*PinvH*H*V;
            Innov=repmat(U'*PinvH*y(:,i),1,L)-Hnq*x;
            
            % Update observation covariance
            R = U' * PinvH * Rfixed * PinvH' * U;
            RinvtInno = R\Innov;
        else % IOPPF ==1, Optimal proposal PF
            %Added by EVV
            Innov1=repmat(y(:,i),1,L)-H*V*x;
            Qnew=V'*Qfixed*V;
            Qpfixed = inv(inv(Qnew)+V'*H'*inv(Rfixed)*H*V);
            QpHRinv = Qpfixed*V'*H'*inv(Rfixed);
            x = x + QpHRinv*Innov1 + mvnrnd(pzeros_physical,Qpfixed,L)';
            
            Innov=U'*PinvH*Innov1;
            % Update observation covariance
            R = U' * PinvH * Rfixed * PinvH' * U;
            Hnq = U'*PinvH*H*V;
            Rnew = R + Hnq*Qnew*Hnq';
            RinvtInno = Rnew\Innov;
            %Rinv = inv(R + Hnq*Qnew*Hnq');
        end
        
        % Reweight
        Tdiag = diag(Innov'*RinvtInno);
        tempering = 1.2; % Tempering usually a little larger than 1.
        Avg=(max(Tdiag)+min(Tdiag))/2;
        Tdiag = (Tdiag-Avg)/tempering;
        LH = exp(-Tdiag/2); %%%% <<<< divided exponent by 2; this is part of the normal distribution
        w=LH.*w;
        %Normalize weights
        [dim,~] = size(w);
        Lones = ones(dim,1);
        w=w/(w'*Lones);
        
        %Resampling (with resamp.m that I provided or using the pseudo code in Peter Jan ,... paper)
        [w,x,NRS,ESS] = resamp(w,x,ResampCutoff);
        Resamps = Resamps + NRS;
        
        if (NRS==1)
            % Projected Resampling
            x = x + V'*(alpha*(U*U') + (1-alpha)*eye(N))*(mvnrnd(zeros(1,N),epsOmega*ones(1,N),L)');
        end
        
    end
    
    % Estimate the truth
    estimate(:,i) = x*w;
    
    if DataProjection ==3 || PhysicalProjection ==3
        if PhysicalProjection == 3
            V0=q;
        end
        %Get AUS projections: [U,V1,V0]=getausproj( ... )
        if PhysicalProjection == 3
            NumLEs=p_physical;
        else
            NumLEs=p_data;
        end
        [q,LE] = getausproj(N,NumLEs,F,t,V*estimate(:,i),h,q,LE);
        if PhysicalProjection == 3
            V=q;
        end
        if DataProjection ==3
            U=q(:,1:numModes_data);
        end
    end
    
    % propogate particles
    x = V'*dp4(F,t,V0*x,h);
    
    % Compare estimate and truth
    diff_orig= truth(:,i) - (V*estimate(:,i));
    diff_proj= (V * V'* truth(:,i)) - (V*estimate(:,i));
    RMSE_orig = sqrt(diff_orig'*diff_orig/N);
    RMSE_proj = sqrt(diff_proj'*diff_proj/N);
    
    RMSEave_orig = RMSEave_orig + RMSE_orig;
    RMSEave_proj = RMSEave_proj + RMSE_proj;
    %Save XC
    xbar=V*estimate(:,i);
    truth_common=truth(:,i);
    Ensbar = mean(xbar);
    Trubar = mean(truth_common);
    XC = (xbar-Ensbar)'*(truth_common-Trubar)/(norm(xbar-Ensbar,2)*norm(truth_common-Trubar));
    
    xbar=V*estimate(:,i);
    truth_common=V*V'*truth(:,i);
    Ensbar = mean(xbar);
    Trubar = mean(truth_common);
    XCproj = (xbar-Ensbar)'*(truth_common-Trubar)/(norm(xbar-Ensbar,2)*norm(truth_common-Trubar));
    
    
    
    % Save to plot
    if mod(i,ObsMult)==0
        %Save RMSE values
        Time(iRMSE)=t;
        RMSEsave(iRMSE)=RMSE_orig;
        RMSEsave_proj(iRMSE)=RMSE_proj;
        XCsave(iRMSE)=XC;
        XCprojsave(iRMSE)=XCproj;
        ESSsave(iRMSE)=ESS;
        iRMSE = iRMSE+1;
    end
    
    t = t+h;
end
%%
RMSEave_orig = RMSEave_orig/Numsteps;
RMSEave_proj = RMSEave_proj/Numsteps;
ResampPercent = ObsMult*Resamps/Numsteps*100;

