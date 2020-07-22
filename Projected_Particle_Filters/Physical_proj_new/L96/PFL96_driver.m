%% Initialization
close all; clear all;clc;
rng(1331);
F = @FLor95; %Physical model
N =100; % N:Original model dimension
% Build Model (via ODE45)
dt=1.E-2; % Model output time step
ModelSteps = 50000; % Number of time steps in building model
T=ModelSteps*dt;
Built_Model= buildModel(N,F,ModelSteps,T);
model_output = Built_Model';
%%
% figure(1)
% % contourf(model_output,'LineStyle','none');
% contourf(model_output(1:40,1:2000),'LineStyle','none');
% colormap(redblue)
% caxis([-0.65, 0.65])
% title('Truth')
% xlabel('Time')
% ylabel('N ')
% % xticklabels(xticks*dt)
%% Type of particle filter
% Use of standard PF or OP-PF (iOPPF=0 => standard PF, iOPPF=1 => OP-PF)
iOPPF=1;

%% Projection_type(0 = no projection, 1 POD, 2 DMD, 3 AUS)
PhysicalProjection =1;
DataProjection = 2;
tolerance_physical = 20; % POD_modes
tolerance_data = 20; % POD_modes
numModes_physical = 20;% DMD_modes, for physical
numModes_data = 20; % DMD_modes, for data

[Ur_physical,p_physical,pzeros_physical] = ...
    Projection_physical_type(PhysicalProjection, numModes_physical,tolerance_physical,N,model_output,dt);
[Ur_data,p_data,pzeros_data] = Projection_data_type(DataProjection ,numModes_data,tolerance_data,N,model_output,dt);

%% Particle Filter Information
L=50;%Number of particles
IC = zeros(N,1);
IC(1)=1; % Particle ICs

alpha=0.35;%alpha value for projected resampling

% Number of computational steps and step size
h=1.E-2;
Numsteps=T/h;

ObsMult=5; % Observe and every ObsMult steps

%Observation Variance
epsR = 0.01;
%Model Variance
epsQ = 0.01;
% IC Variance
epsOmega =0.0027;
%Initial condition
epsIC = 0.01;
%Observe every inth variable.
inth=1;
%Call Init
[M,H,PinvH,IC,q,LE,w,R,Rinv,Q,Omega,ICcov,Lones,Mzeros] = ...
    Init_L96(F,IC,h,N,inth,ModelSteps,p_physical,L,epsR,epsQ,epsOmega,epsIC);

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
iRMSE=1;
[M,H,PinvH] = new_Init_L96(N,inth,V);
Q=V'*Q*V;
x=V'*u;
%%save to plot the RMSE on original and projected space
% diff_plot=[];
% diff_proj_plot=[];
for i=1:Numsteps
    % Estimate the truth    
    estimate(:,i) = x*w;
    
    % Get projection of the Data Model
    [U] = projectionToggle_data(DataProjection,N,Ur_data,p_data); 
    % Update noise covariance
    R = U' * PinvH * Rfixed * PinvH' * U;
    Rinv = inv(R);
    
    if mod(i,ObsMult)==0
        %At observation times, Update weights via likelihood, add noise
        if (iOPPF==0) % Standard Particle Filter             
            x = x + mvnrnd(pzeros_physical,Q,L)';            
            Hnq = U'*PinvH*H*V;
            Innov=repmat(U'*PinvH*y(:,i),1,L)-Hnq*x;
            
        else % IOPPF ==1, Optimal proposal PF
            Hnq = U'*PinvH*H*V; 
            Qpinv = inv(Q) + Hnq'*Rinv*Hnq;
            Qp = inv(Qpinv);
            Innov=repmat(U'*PinvH*y(:,i),1,L)-Hnq*x;
            x = x + Qp*Hnq'*Rinv*Innov + mvnrnd(pzeros_physical,Qp,L)';
            Rinv = inv(R + Hnq*Q*Hnq');
        end
                
        % Reweight
        Tdiag = diag(Innov'*Rinv*Innov);
        tempering = 1.2; % Tempering usually a little larger than 1.
        Avg=(max(Tdiag)+min(Tdiag))/2;
        Tdiag = (Tdiag-Avg)/tempering;
        
        Tdiag = -Tdiag/2;
        logw = Tdiag + log(w);
               
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
        
        
        %Resampling (with resamp.m that I provided or using the pseudo code in Peter Jan ,... paper)
        [w,x,NRS] = resamp(w,x,0.5);
        Resamps = Resamps + NRS;
               
        if (NRS==1)
            % Projected Resampling
            x = x + V'*(alpha*(U*U') + (1-alpha)*eye(N,1))*(mvnrnd(zeros(N,1),Omega,L)');
        end
        
    end
    
    % propogate particles
    x = V'*dp4(F,t,V*x,h);

    % Compare estimate and truth 
    diff_orig= truth(:,i) - (V*estimate(:,i));
    diff_proj= (V * V'* truth(:,i)) - (V*estimate(:,i));
%     diff_plot(:,i)= truth(:,i) - (V*estimate(:,i)); %%to plot the difference on original space 
%     diff_proj_plot(:,i)= (V * V'* truth(:,i)) - (V*estimate(:,i));%%to plot the difference on projected space 
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


% figure(2)
% contourf(diff_plot,'LineStyle','none')
% colormap(redblue)
% caxis([-0.65, 0.65])
% title('The difference in the Original Space')
% xlabel('Time')
% ylabel('N ')

% figure(3)
% contourf(diff_proj_plot,'LineStyle','none')
% colormap(redblue)
% colorbar;
% caxis([-0.65, 0.65])
% title('The difference in the Projected Space')
% xlabel('Time')
% ylabel('N ')

%%
% % epsRR=epsR*ones(1,length(RMSEsave));
% % loyolagreen = 1/255*[0,104,87];
% figure(4)
% plot(Time,RMSEsave, 'b--', 'LineWidth', 1.5)
% grid on
% hold on;
% plot(Time,RMSEsave_proj,'r--','LineWidth', 1.5)
% plot(Time,epsRR,':','Color', loyolagreen,'LineWidth', 1.5)
% % plot(Time,epsRR,'m:','LineWidth', 1.5)
% xlabel('Time')
% ylabel('RMSE')
% ylim([0 0.3])
% % title('POD Projection')
% % title('DMD Projection')
% % title('POD and DMD Projection')
% % legend('RMSE','Observation error','Location', 'Best')
% legend('RMSE Original','RMSE Projected','Observation error','Location', 'Best')


%%
RMSEave_orig = RMSEave_orig/Numsteps
RMSEave_proj = RMSEave_proj/Numsteps
ResampPercent = ObsMult*Resamps/Numsteps
