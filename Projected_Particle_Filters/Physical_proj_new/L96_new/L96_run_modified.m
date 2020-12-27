close all; clear;clc;
ObsMult=5; % % Observe and every ObsMult steps(10 with F=3.5, 5 with F=8)
epsQ=1E-1;
epsR=1E-2;
N=400;
Num=40;
Mult=N/Num;
%% Projection_type(0 = no projection, 1 POD, 2 DMD, 3 AUS)
PhysicalProjection =2;
DataProjection =2 ;
Num_trials=10;
for j=1:Num
    for k=1:Num_trials
        numModes_physical=j*Mult+1;%DMD
        tolerance_physical=j*Mult;%POD
        nummodes=numModes_physical-1;
        if nummodes==N
            PhysicalProjection=0;
        end
        [Time(:,j,k),RMSEsave(:,j,k), RMSEsave_proj(:,j,k), XCsave(:,j,k), XCprojsave(:,j,k), ESSsave(:,j,k), ResampPercent(:,j,k)]= L96_modified...
            (numModes_physical,epsQ, epsR,tolerance_physical,PhysicalProjection,DataProjection,ObsMult,N);
    end
end
%% Save to mat file
PhysicalProjection=1;
filename = sprintf('L968_%2d%2d_%2d_%2d_%4d.mat',PhysicalProjection,DataProjection,epsQ,epsR,N)
params.PhysicalProjection = PhysicalProjection;
params.DataProjection = DataProjection;
params.epsQ = epsQ;
params.epsR=epsR;
params.numModes=numModes_physical;
params.Num=Num;
params.N=N;
params.Mult=Mult;

results.RMSEsave = RMSEsave;
results.RMSEsave_proj = RMSEsave_proj;
results.XCsave = XCsave;
results.XCsave_proj =XCprojsave;
results.ESSsave= ESSsave;
results.ResampPercent =ResampPercent;
results.Time =Time;

save(filename,'params','results');
