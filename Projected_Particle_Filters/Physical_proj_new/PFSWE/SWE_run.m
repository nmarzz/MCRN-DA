close all; clear;clc;

%% SET ALL PARAMETERS HERE
epsOmega =0.0000001; %For inth = 1
%Observe every inth variable.
Numsteps=1441;
inth=100;%observations run for 10,100,1000 
epsQ=1; %Run for epsQ=1
epsR=1E-2;
Num=1; % number of model reduction orders
L=5; % numbe3r of particles

% Observed variables scenario following Paulina et al.
% scenario 1: observe inth u and v's
% scenario 2: observe inth everything
% scenario 3: observe inth h
scenario =2 ;
%% Type of particle filter
% Use of standard PF or OP-PF (iOPPF=0 => standard PF, iOPPF=1 => OP-PF)
iOPPF=1;

%% Projection_type(0 = no projection, 1 POD, 2 DMD, 3 AUS)
PhysicalProjection =0;
DataProjection =0;

Num_trials=10;

Mult=100/Num; % don't touch!


%% Parameter storage
filename = sprintf('SWE_%s_%s_Q%1.2f_R%1.2f_inth%04d_S%d_L%d.mat',...
    proj2string(PhysicalProjection),...
    proj2string(DataProjection),...
    epsQ,...
    epsR,...
    inth,...
    scenario,...
    L)

params.PhysicalProjection = PhysicalProjection;
params.DataProjection = DataProjection;
params.epsQ = epsQ;
params.epsR=epsR;
params.iOPPF = iOPPF;
params.Num=Num;
params.Mult=Mult;
params.epsOmega=epsOmega;
params.Num_trials = Num_trials;

%% Main computation loop
disp('SWE Main loop started');

for j=1:Num
    numModes_physical=j*Mult+1;%DMD
    tolerance_physical=j*Mult;%POD
    nummodes=numModes_physical-1;
    params.numModes=numModes_physical; % update to last simulation that was reached for storage purposes
    
    for k=1:Num_trials % REPLACE WITH PLAIN FOR IF YOU DON'T WANT PARALLEL EXECUTION
        [Time(:,j,k),RMSEsave(:,j,k),ResampPercent(:,j,k),ESSsave(:,j,k)]=PFSWErun...
            (numModes_physical,epsQ,epsR,tolerance_physical,iOPPF,PhysicalProjection,DataProjection,scenario,epsOmega, inth,Numsteps,L);
    end
    disp("SWE CHECKPOINT: Mode choice " +  j + "/" + Num + " - Trials: " + Num_trials + " - TOTAL RUNTIME " + string( duration([0, 0, toc]) )); 

end

%% Results storage

results.RMSEsave = RMSEsave;
% results.RMSEsave_proj = RMSEsave_proj;
results.ResampPercent =ResampPercent;
results.ESS=ESSsave;
results.Time =Time;

%% Save to drive
save(filename,'params','results');

%% AUXILIARY FUNCTION
function out = proj2string(val)
switch(val)
    case 0
        out = 'NON';
    case 1
        out = 'POD';
    case 2
        out = 'DMD';
    case 3
        out = 'AUS';
    otherwise
        error("Unknown method")
end
end


