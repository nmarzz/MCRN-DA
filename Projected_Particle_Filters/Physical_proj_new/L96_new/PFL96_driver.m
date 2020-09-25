%% Initialization
close all; clear all;clc;
%rng(1331);
rng(1330);
F = @FLor95; %Physical model
%N =100; % N:Original model dimension
N =40; % N:Original model dimension
% Build Model (via ODE45)
dt=1.E-2; % Model output time step
%ModelSteps = 50000; % Number of time steps in building model
ModelSteps = 100000; % Number of time steps in building model
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
PhysicalProjection =3;
DataProjection = 3;
tolerance_physical = 10; % POD_modes
tolerance_data = 5; % POD_modes
%numModes_physical = 20;% DMD_modes/AUS_modes, for physical
%numModes_data = 20; % DMD_modes/AUS_modes, for data
numModes_physical = 8;% DMD_modes/AUS_modes, for physical
numModes_data = 4; % DMD_modes/AUS_modes, for data

[Ur_physical,p_physical,pzeros_physical] = ...
    Projection_physical_type(PhysicalProjection, numModes_physical,tolerance_physical,N,model_output,dt);
[Ur_data,p_data,pzeros_data] = Projection_data_type(DataProjection ,numModes_data,tolerance_data,N,model_output,dt);

%% Particle Filter Information
L=20;%Number of particles
IC = zeros(N,1);
IC(1)=1; % Particle ICs

%alpha=0.35;%alpha value for projected resampling
alpha=0.99;%alpha value for projected resampling
%ResampCutoff=0.3;
ResampCutoff=0.5;

% Number of computational steps and step size
h=1.E-2;
Numsteps=T/h;

ObsMult=5; % Observe and every ObsMult steps

%Observation Variance
%epsR = 0.1;
epsR = 0.01;
%epsR = 0.25;
%Model Variance
%epsQ = 0.1;
epsQ = 1;
%epsQ = 0.01;
% IC Variance
%epsOmega =0.0027;
epsOmega =0.0001;
%epsOmega =0.01;
%epsOmega =0.00001;
%Initial condition
%epsIC = 0.01;
epsIC = epsQ;
%Observe every inth variable.
%inth=2;
inth=1;
%Call Init
p_physical
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
iRMSE=1;
[M,H,PinvH] = new_Init_L96(N,inth,V);
x=V'*u;

Qfixed=Q;
Q=V'*Q*V;

%Added by EVV
Qnew=V'*Qfixed*V;
Qpfixed = inv(inv(Qnew)+V'*H'*inv(Rfixed)*H*V);
QpHRinv = Qpfixed*V'*H'*inv(Rfixed);

%%save to plot the RMSE on original and projected space
% diff_plot=[];
% diff_proj_plot=[];
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
            %Rinv = inv(R);
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

        %NEW: start
        LH = exp(-Tdiag/2); %%%% <<<< divided exponent by 2; this is part of the normal distribution
        w=LH.*w;

        %Normalize weights
        [dim,~] = size(w);
        Lones = ones(dim,1);
        w=w/(w'*Lones);
        %NEW: end
        
        %Resampling (with resamp.m that I provided or using the pseudo code in Peter Jan ,... paper)
        [w,x,NRS,ESS] = resamp(w,x,ResampCutoff);
        Resamps = Resamps + NRS;
               
        if (NRS==1)
            % Projected Resampling
            x = x + V'*(alpha*(U*U') + (1-alpha)*eye(N))*(mvnrnd(zeros(1,N),epsOmega*ones(1,N),L)');
            %x = x + V'*(alpha*(U*U') + (1-alpha)*eye(N,1))*(mvnrnd(zeros(N,1),Omega,L)');
            %x = x + (alpha*(U*U') + (1-alpha)*eye(N,1))*(mvnrnd(zeros(N,1),Omega,L)');
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
    %x = dp4(F,t,x,h);

    % Compare estimate and truth 
    diff_orig= truth(:,i) - (V*estimate(:,i));
    diff_proj= (V * V'* truth(:,i)) - (V*estimate(:,i));
    %diff_orig= truth(:,i) - (estimate(:,i));
    %diff_proj= (V * V')* diff_orig;
%     diff_plot(:,i)= truth(:,i) - (V*estimate(:,i)); %%to plot the difference on original space 
%     diff_proj_plot(:,i)= (V * V'* truth(:,i)) - (V*estimate(:,i));%%to plot the difference on projected space 
    RMSE_orig = sqrt(diff_orig'*diff_orig/N);
    RMSE_proj = sqrt(diff_proj'*diff_proj/N);
    MAE_orig = (sum(abs(diff_orig)))/N;
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

%Display the RMSE vs. time
[~,Steps] = size(Time)
ObsErr = linspace(sqrt(epsR),sqrt(epsR),Steps);
ModErr = linspace(sqrt(epsQ),sqrt(epsQ),Steps);
figure(2)
plot(Time,RMSEsave)
hold on
plot(Time,RMSEsave_proj)
plot(Time, ObsErr,'k--')
plot(Time,ModErr,'r--')
%axis([t0 tf 0 2])
legend('RMSE','RMSE projected','Observation Error','Model Error','Location','NorthWest')
xlabel('Time')
ylabel('RMSE')
title('Root Mean Squared Error')

figure(3)
XCave=mean(XCsave)
plot(Time,XCsave)
hold on
plot(Time, XCprojsave,'k--')
legend('XC','XCproj')
xlabel('Time')
ylabel('XC')
title('Pattern Correlation Coefficient')


figure(4)
ESSave=mean(ESSsave)
plot(Time,ESSsave)
hold on
Cutoff = linspace(L*ResampCutoff,L*ResampCutoff,Steps);
plot(Time, Cutoff,'k--')
xlabel('Time')
ylabel('ESS')
title('Effective Sample Size')

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
ResampPercent = ObsMult*Resamps/Numsteps*100
