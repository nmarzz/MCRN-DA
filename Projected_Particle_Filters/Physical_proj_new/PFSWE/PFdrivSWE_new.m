%% Initialization
clear all;clc;close all;
%rng(1331);
rng(1330);
% SWE preamble
%load('SWE_run_1day.mat');
% load('SWE_run_2days.mat');
% load('SWE_run_3days.mat');
load('SWE_run_4days.mat');
parsanim.H = zeros(pars.nx, pars.ny);
parsanim.x=(0:pars.nx-1).*pars.dx; % Zonal distance coordinate (m)
parsanim.y=(0:pars.ny-1).*pars.dy; % Meridional distance coordinate (m)
F = @(t,x) formod(t,x,dt,pars);
Built_Model= x_save;
N =length(Built_Model);
IC = Built_Model(:,(end-1)/2);
% figure(5)
% contourf(Built_Model(15000:15500,1000:1500),'LineStyle','none');
% colormap(redblue)
% % cmocean('balance')
% caxis([-6, 6])
% title('Truth')
% xlabel('Time')
% ylabel('N ')
%% Type of particle filter
% Use of standard PF or OP-PF (iOPPF=0 => standard PF, iOPPF=1 => OP-PF)
iOPPF=1;

%% Projection_type(0 = no projection, 1 POD, 2 DMD, 3 AUS)
PhysicalProjection = 2;
DataProjection = 2;
tolerance_physical =20; % POD_modes
tolerance_data = 10; % POD_modes
numModes_physical = 30;% DMD_modes, for physical
numModes_data = 20; % DMD_modes, for data

%model_output = Built_Model;
model_output = Built_Model(:,(end-1)/4:(end-1)*3/4);
[Ur_physical,p_physical,pzeros_physical] = ...
    Projection_physical_type(PhysicalProjection ,numModes_physical,tolerance_physical,N,model_output,dt);
[Ur_data,p_data,pzeros_data] = Projection_data_type(DataProjection ,numModes_data,tolerance_data,N,model_output,dt);

%% Particle Filter Information
L=20;%Number of particles
% L=20;%Number of particles
%alpha =0;%alpha value for projected resampling
alpha =.99;%alpha value for projected resampling
ResampCutoff = 0.3;
% Number of computational steps and step size
ObsMult=1; % Observe and every ObsMult steps
h = dt/ObsMult;
%Numsteps = size(Built_Model,2)*ObsMult;
%Numsteps =Numsteps/2;
Numsteps=1000;
NumstepsBig=size(Built_Model,2);

%Observation Variance
%alph = 0.01;
alph = 1;
bet = 0.01;
%bet = .1;
%epsR = 0.01;
epsR = bet;
%Model Variance
%epsQ = 0.1;
epsQ = alph;
% IC Variance
%epsOmega =0.0027;
epsOmega =0.0000001; %For inth = 1
% epsOmega =0.001; %For inth = 1000
%Initial condition
%epsIC = 0.01;
epsIC = 1;
%Observe every inth variable.
inth=1;
% inth=1000;
%Call Init
[M,IC,wt,R,Rinv,Q,Omega,ICcov,Lones,Mzeros] = Init_simp(IC,N,inth,L,epsR,epsQ,epsOmega,epsIC,alph, bet);
%Init(IC,N,inth,L,epsR,epsQ,epsOmega,epsIC);

%Add noise N(0,ICcov) to ICs to form different particles

%u = repmat(IC,1,L) + normrnd(0,ICcov,N,L); % Noise for IC
u = repmat(IC,1,L) + mvnrnd(zeros(1,N),ICcov*ones(1,N),L)'; % Noise for IC

%% Generate observations from "Truth"

y=zeros(size(IC(1:inth:end),1),NumstepsBig);
gen_ics = IC;
t = 0;
for i = 1:NumstepsBig
    truth(:,i) = gen_ics;
    gen_ics = formod(t,gen_ics,dt,pars);
    if mod(i,ObsMult) == 0
        y(:,i) = truth(1:inth:end,i);
    end
    t = t + h;
end
%y = y + normrnd(0,R,M,size(Built_Model,2));
%y = y + mvnrnd(Mzeros',R*ones(1,M),Numsteps)'; % Add noise to observations
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
RMSEave_proj=0;
RMSEave_relRMSE=0;
iRMSE=1;

Q=V'*Q*V; %with projection
Qnew=Q;

x=V'*u;

%Added by EVV
HV=Hx(V,inth);  % with our assumptions left multipling
UPinvH=Hx(U,inth)';

Qpfixed = inv(inv(Qnew)+HV'*inv(Rfixed)*HV);
QpHRinv = Qpfixed*HV'*inv(Rfixed);

%p_physical = 10;
% diff_plot=[];
% diff_proj_plot=[];
for i=1:Numsteps
    t
    %[U] = projectionToggle_data(DataProjection,Ur_data,p_data);
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
                x = x + QpHRinv*Innov1 + mvnrnd(pzeros_physical',Qpfixed*ones(1,sizepzp),L)';
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
        
        %NEW: start
        LH = exp(-Tdiag/2); %%%% <<<< divided exponent by 2; this is part of the normal distribution
        wt=LH.*wt;
        
        %Normalize weights
        [dim,~] = size(wt);
        Lones = ones(dim,1);
        wt=wt/(wt'*Lones);
        %NEW: end
        
        
        %        logw = Tdiag + log(wt);
        %
        %        [~,idx] = min(abs(logw-((max(logw) - min(logw))/2))); % find index of weight closest to middle value
        %        logw([1 idx]) = logw([idx 1]);
        %        x(:,[1 idx]) = x(:,[idx 1]);
        %        toEXP = logw(2:end) - logw(1);
        %        toEXP=min(toEXP,709);
        %        toEXP=max(toEXP,-709);
        %        toSum = exp(toEXP);
        %        normalizer = logw(1) + log1p(sum(toSum));
        %        logw = logw - normalizer;
        %        wt = exp(logw);
        
        %Resampling (with resamp.m that I provided or using the pseudo code in Peter Jan ,... paper)
        [wt,x,NRS] = resamp(wt,x,ResampCutoff);
        Resamps = Resamps + NRS;
        if (NRS==1)
            % Projected Resampling
            x = x + (alpha*(V'*U*U') + (1-alpha)*V')*(mvnrnd(zeros(1,N),epsOmega*ones(1,N),L)');
            %x = x + V'*(alpha*(U*U') + (1-alpha))*(mvnrnd(zeros(1,N),epsOmega*ones(1,N),L)');
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
    diff_proj= V*(V'* truth(:,i) - estimate(:,i));
    %     diff_plot(:,i)= truth(:,i) - (V*estimate(:,i));
    %diff_proj_plot(:,i)= V*(V'* truth(:,i) - estimate(:,i));
    RMSE_orig = sqrt(diff_orig'*diff_orig/N)
    RMSE_proj = sqrt(diff_proj'*diff_proj/N)
    %relRMSE_orig=sqrt(diff_orig'*diff_orig/N)/(sqrt(((truth(:,i))'*truth(:,i))/N));
    %MAE_orig = (sum(abs(diff_orig)))/N;
    RMSEave_orig = RMSEave_orig + RMSE_orig;
    RMSEave_proj = RMSEave_proj + RMSE_proj;
    %RMSEave_relRMSE =RMSEave_relRMSE+ relRMSE_orig;
    
%     figure(1)
%     animate_NEW2(t,V*estimate(:,i),pars,parsanim,'Mean')
%     
%     figure(2)
%     animate_NEW2(t,truth(:,i),pars,parsanim,'Truth')
%     
%     figure(3)
%     animate_NEW2(t,diff_orig,pars,parsanim,'Difference (original)')
%     
%     figure(4)
%     animate_NEW2(t,diff_proj,pars,parsanim,'Difference (projected)')
    
    % Save to plot
    if mod(i,ObsMult)==0
        %Save RMSE values
        Time(iRMSE)=t;
        RMSEsave(iRMSE)=RMSE_orig;
        RMSEsave_proj(iRMSE)=RMSE_proj;
        %RMSEsave_relRMSE(iRMSE) = relRMSE_orig;
        
        iRMSE = iRMSE+1;
    end
    
    t = t+h;
end
%
% % epsRR=epsR*ones(1,length(RMSEsave));
% % % loyolagreen = 1/255*[0,104,87];
figure
plot(Time,RMSEsave, 'r-.', Time,RMSEsave_proj, 'g-.',Time,sqrt(bet)*ones(size(Time,2),1), 'k-.','LineWidth', 2)
grid on
% hold on;
% plot(Time,RMSEsave_proj,'b-.','LineWidth', 1.5)
% plot(Time,epsRR,'g-.','LineWidth', 1.5)
xlabel('Time')
ylabel('RMSE')
legend('Original','Projected','Obs Error')
% ylim([0 60])
% title('POD Projection(r=15)')
% title('Stable RMSE')
% legend('RMSE Original','RMSE Projected','Location', 'Best')
% % % %%
% figure(2)
% contourf(diff_plot(15000:15500,1000:1500),'LineStyle','none')
% colormap(redblue)
% % % colorbar
% caxis([-6, 6])
% title('The difference in the Original Space')
% xlabel('Time')
% ylabel('N ')
% %
% % figure(3)
% % contourf(diff_proj_plot(15000:15500,1000:1500),'LineStyle','none')
% % colormap(redblue)
% % colorbar;
% % caxis([-6, 6])
% % title('The difference in the Projected Space')
% % xlabel('Time')
% % % ylabel('N ')
%
RMSEave_orig = RMSEave_orig/Numsteps
RMSEave_proj = RMSEave_proj/Numsteps
% RMSEave_relRMSE=RMSEave_relRMSE/Numsteps;
ResampPercent = ObsMult*Resamps/Numsteps*100
