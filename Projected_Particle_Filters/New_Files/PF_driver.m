%% Initialization
clear all;clc;
%Projection Parameters
%(0 = no projection, 1 POD, 2 DMD, 3 AUS)
%Build Model (dimension, model)
Model_Dimension = 40;
Fmod = @FLor95;
Built_Model = buildModel(Model_Dimension,@FLor95);
PhysicalProjection =0;
DataProjection = 1;
Ur_physical=0;
Ur_data=0;
if PhysicalProjection ==1
    %POD
    tolerance = 0.0001;
    Ur_physical= buildPOD(tolerance,Built_Model);
elseif PhysicalProjection ==2
    %DMD
    numModes=50;%number of DMD_modes you want to use
    Ur_physical=buildDMD(numModes,Built_Model);
end
%%
if DataProjection ==1
    %POD
    tolerance = 0.0001;
    Ur_data= buildPOD(tolerance,Built_Model);
elseif DataProjection ==2
    %DMD
    numModes=50;%number of DMD_modes you want to use
    Ur_data=buildDMD(numModes,Built_Model);
end
%% Particle Filter Information
% Type of particle filter
%Use of standard PF or OP-PF (iOPPF=0 => standard PF, iOPPF=1 => OP-PF)
iOPPF=1;
%Number of particles
L=50;
%alpha value for projected resampling
alpha=1;
%Number of time steps
Numsteps = 1000;
%Multiple of the step size for observation time
ObsMult=10;
%Rank of projection, number of Lyapunov exponents for AUS projection
p=10;
%ICs for particles
IC = zeros(Model_Dimension,1); % N -> Model_Dimension incase this causes problems
IC(1)=1;
%Computational time step
h=5.E-3;
%For diagonal (alpha*I) covariance matrices.
%Observation
epsR = 0.01;
%Model
epsSig = 0.01;
%Resampling
epsOmega = 2.E-3;
%Initial condition
epsIC = 0.1;
%Observe every inth variable.
inth=2;
%Call Init
[M,H,PinvH,IC,q,LE,w,R,Rinv,Sig,Omega,ICcov,Lones,Mzeros,Nzeros] = ...
    Init(Fmod,IC,h,Model_Dimension,inth,Numsteps,p,L,epsR,epsSig,epsOmega,epsIC);
Rinvfixed=Rinv;
%Add noise N(0,ICcov) to ICs to form different particles
x = zeros(Model_Dimension,L);
x = repmat(IC,1,L) + mvnrnd(Nzeros,ICcov,L)'; %ICchol*randn(N,L);
estimate(:,1) = x*w;

y=zeros(M,Numsteps);

t=0;
%Generate observations from "Truth"
for i = 1:Numsteps
    truth(:,i) = IC;
    if mod(i,ObsMult)==0
        y(:,i)=H*IC + mvnrnd(Mzeros,R,1)'; %Rchol*rand(M,1); % + Noise from N(0,R)
    end
    IC = dp4(Fmod,t,IC,h);
    t = t+h;
end

%Initial time and time step
t=0;
t0=t;
Resamps=0;
RMSEave=0;
iRMSE=1;
[q_physical] =projectionToggle_Physical(PhysicalProjection,Model_Dimension,Ur_physical,p);
%Loop over observation times
[M,H,PinvH] = new_Init(Model_Dimension,inth,q_physical);
Sig=q_physical'*Sig*q_physical;
x=q_physical'*x;
for i=1:Numsteps
    % Perform projection of the Data Model
    [q_data] = projectionToggle_data(DataProjection,Model_Dimension,Ur_data,p); %chooses which q projection we want
    est=estimate(:,i);
    if mod(i,ObsMult)==0
        %At observation times, Update weights via likelihood
        if (iOPPF==0)
            %Add noise only at observation times
            x = x + mvnrnd(Nzeros,Sig,L)';
            %Standard Particle Filter
            %Standard PF (no projection)
            Innov = repmat(y(:,i),1,L) - H*x;
        else % Ioppf ==1
            Qpinv = inv(Sig) + H'*Rinvfixed*H;
            Qp = inv(Qpinv);
            %Optimal proposal PF
            Innov = repmat(y(:,i),1,L) - H*x;
            x = x + Qp*H'*Rinv*Innov + mvnrnd(Nzeros,Qp,L)';
            Rinv = inv(R + H*Sig*H');
        end
        
        Tdiag = diag(Innov'*Rinv*Innov);
        tempering = 1.2; %%%% <<< including new parameter here for visibility. Tempering usually a little larger than 1.
        Tdiag = (Tdiag-max(Tdiag))/tempering; %%%%% <<<< Think dividing the exponent is dangerous; this was tempering with an unknown coefficient.
        LH = exp(-Tdiag/2); %%%% <<<< divided exponent by 2; this is part of the normal distribution
        w=LH.*w;
        %Normalize weights
        w=w/(w'*Lones);
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
    diff = truth(:,i)-estimate(:,i);
    RMSE = sqrt(diff'*diff/Model_Dimension)
    RMSEave = RMSEave + RMSE;
    
    if mod(i,ObsMult)==0
        %Save RMSE values
        Time(iRMSE)=t;
        RMSEsave(iRMSE)=RMSE;
        iRMSE = iRMSE+1;
        
        %Plot
        % yvars=colon(1,inth,N);
        % vars = linspace(1,N,N);
        % sz=zeros(N,1);
        % plots(1) = plot(vars,truth(:,i),'ro-');
        % hold on
        % plots(2) = plot(vars,estimate(:,i+1),'bo-');
        % for j=1:L
        %   sz(:)=w(j)*80*L;
        %   scatter(vars,x(:,j),sz,'b','filled');
        % end
        % plots(3) = plot(yvars,y(:,i),'g*','MarkerSize',20);
        % title(['Time = ',num2str(t)])
        % legend(plots(1:3),'Truth','Estimate','Obs');
        % pause(1);
        % hold off
        %
    end
    t = t+h;
end
figure(2)
plot(Time,RMSEsave);
RMSEave = RMSEave/Numsteps
ResampPercent = ObsMult*Resamps/Numsteps

LE = LE/(t-t0)