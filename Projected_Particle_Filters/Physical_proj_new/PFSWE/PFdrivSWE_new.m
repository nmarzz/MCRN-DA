%% Initialization
clear all;clc;
rng(1330);
% SWE preamble
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
iOPPF=0;
%% Projection_type(0 = no projection, 1 POD, 2 DMD, 3 AUS)
PhysicalProjection = 2;
DataProjection =1;
tolerance_physical =130; % POD_modes
tolerance_data = 10; % POD_modes
numModes_physical =130;% DMD_modes, for physical
numModes_data =10; % DMD_modes, for data

%model_output = Built_Model;
model_output = Built_Model(:,(end-1)/4:(end-1)*3/4);
[Ur_physical,p_physical,pzeros_physical] = ...
    Projection_physical_type(PhysicalProjection ,numModes_physical,tolerance_physical,N,model_output,dt);
[Ur_data,p_data,pzeros_data] = Projection_data_type(DataProjection ,numModes_data,tolerance_data,N,model_output,dt);

%% Particle Filter Information
L=20;%Number of particles
ResampCutoff = 0.3;
% Number of computational steps and step size
ObsMult=1; % Observe and every ObsMult steps
h = dt/ObsMult;
Numsteps=500;
NumstepsBig=size(Built_Model,2);
% %% FOR PF
%Observation Variance
alpha =0.99;%alpha value for projected resampling
alph = 0.001;%PF
bet =1;
epsR = bet;
%Model Variance
epsQ = alph;
%Initial condition
epsIC =0.01 ;


% %% For PF-OP
% alpha =.99;%alpha value for projected resampling
% alph = 1;
% bet = 0.01;
% epsR = bet;
% %Model Variance
% epsQ = alph;
% %Initial condition
% epsIC = 0.01;
%%
% IC Variance
epsOmega =0.0000001; %For inth = 1
% epsOmega =0.001; %For inth = 1000
%Observe every inth variable.
% inth=1000;
inth=1;
%Call Init
[M,IC,wt,R,Rinv,Q,Omega,ICcov,Lones,Mzeros] = Init_simp(IC,N,inth,L,epsR,epsQ,epsOmega,epsIC,alph, bet);
%Add noise N(0,ICcov) to ICs to form different particles
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
% RMSEave_relRMSE=0;
iRMSE=1;
XC_save_ave=0;
Q=V'*Q*V; %with projection
Qnew=Q;

x=V'*u;

%Added by EVV
HV=Hx(V,inth);  % with our assumptions left multipling
UPinvH=Hx(U,inth)';

Qpfixed = inv(inv(Qnew)+HV'*inv(Rfixed)*HV);
QpHRinv = Qpfixed*HV'*inv(Rfixed);
% diff_plot=[];
% diff_proj_plot=[];
% XC_save=[];
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

        
        %Take log of the updated weights
        logw = Tdiag + log(wt);
        %To avoid underflow and overflow
        logw=min(logw,709);
        logw=max(logw,-709);
        %Exponentiate and normalize
        wt = exp(logw);
        wt = wt/sum(wt);
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
    RMSE_orig = sqrt(diff_orig'*diff_orig/N)
    RMSE_proj = sqrt(diff_proj'*diff_proj/N)
    Ensbar = mean(V*estimate(:,i));
    Trubar = mean( truth(:,i));
    XC_save= (V*estimate(:,i)-Ensbar)'*(truth(:,i)-Trubar)/(norm(V*estimate(:,i)-Ensbar,2)*norm( truth(:,i)-Trubar,2));
    %     diff_plot(:,i)= truth(:,i) - (V*estimate(:,i));
    %diff_proj_plot(:,i)= V*(V'* truth(:,i) - estimate(:,i));
    
    %relRMSE_orig=sqrt(diff_orig'*diff_orig/N)/(sqrt(((truth(:,i))'*truth(:,i))/N));
    %MAE_orig = (sum(abs(diff_orig)))/N;
    
    RMSEave_orig = RMSEave_orig + RMSE_orig;
    RMSEave_proj = RMSEave_proj + RMSE_proj;
%     XC_save_2= XC_save_2+ XC_save;
    %RMSEave_relRMSE =RMSEave_relRMSE+ relRMSE_orig;
    
    %     figure(1)
    %     animate(t,V*estimate(:,i),pars,parsanim,'Mean')
    %
    %     figure(2)
    %     animate(t,truth(:,i),pars,parsanim,'Truth')
    %
    %     figure(3)
    %     animate(t,diff_orig,pars,parsanim,'Difference (original)')
    %
    %     figure(4)
    %     animate(t,diff_proj,pars,parsanim,'Difference (projected)')
    
    % Save to plot
    if mod(i,ObsMult)==0
        %Save RMSE values
        Time(iRMSE)=t;
        RMSEsave(iRMSE)=RMSE_orig;
        RMSEsave_proj(iRMSE)=RMSE_proj;
        XC_save_ave(iRMSE)=XC_save;
        iRMSE = iRMSE+1;
    end
    
    t = t+h;
end
% % No_proj=load('No_proj_RMSE.mat');%PF-OP
No_proj=load('No_proj_PF.mat');%PF
RMSE_no_proj=No_proj.RMSEsave;
TOLC=ptc12(9);
figure(1)
semilogy(Time,RMSEsave,'Color', TOLC(1,:),'LineStyle','-','LineWidth', 2)
hold on
semilogy(Time,RMSEsave_proj,'Color', TOLC(2,:),'LineStyle','--','Marker','+','LineWidth', 2)
semilogy(Time,RMSE_no_proj,'Color', TOLC(7,:),'LineStyle',':','Marker','.','LineWidth',2)
semilogy(Time,sqrt(bet)*ones(size(Time,2),1),'k-.','LineWidth', 2)
grid on
xlabel('Time','fontsize',14,'interpreter','latex','FontName', 'Times New Roman','fontweight','bold')
ylabel('RMSE','fontsize',14,'interpreter','latex','FontName', 'Times New Roman','fontweight','bold')
legend('Model Space','Projected Space','No Reduction','Observation Error','Location', 'Best','fontsize',13,'interpreter','latex','FontName', 'Times New Roman','fontweight','bold')
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
 hold off
% figure(2)
% plot(Time,real(XC_save_ave),'Color', TOLC(1,:),'LineStyle','-','LineWidth', 2)
% grid on
% xlabel('Time','fontsize',14,'interpreter','latex','FontName', 'Times New Roman','fontweight','bold')
% ylabel('Pattern Correlations ','fontsize',14,'interpreter','latex','FontName', 'Times New Roman','fontweight','bold')
% % yticks ([0.9 0.99 0.9999 1])
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)

RMSEave_orig = RMSEave_orig/Numsteps
RMSEave_proj = RMSEave_proj/Numsteps
% RMSEave_relRMSE=RMSEave_relRMSE/Numsteps;
ResampPercent = ObsMult*Resamps/Numsteps*100
