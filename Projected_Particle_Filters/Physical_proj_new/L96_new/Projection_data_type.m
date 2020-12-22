function [Ur_data,p,Nzeros] = Projection_data_type(DataProjection,numModes,tolerance,Model_Dimension,Built_Model,dt )

%Physical Projection
if DataProjection == 0
    Ur_data=0;
    p=Model_Dimension;
    Nzeros=zeros(Model_Dimension,1);
elseif DataProjection == 1
    %POD
    Ur_data= buildPOD(tolerance,Built_Model);
    p = size(Ur_data,2);
    Nzeros=zeros(p,1);
elseif DataProjection == 2
    %DMD
    DataMatrix=Built_Model;%snapshot matrix
    Nsteps = size(DataMatrix, 2); % number of columns
    stepVALUE = round( linspace(1, Nsteps, 5) ); % steps at which b coefficients will be computed
    VALUEmean=true;
    VALUEsortbyb=true;
    M = mean(DataMatrix, 2); % 2 = mean across columns
    r=numModes;
    out1 = dmd( DataMatrix, dt, r-1, 'removemean', VALUEmean,'step',stepVALUE,'sortbyb', VALUEsortbyb);
    bM = norm(M);
    Mnormalized = M/bM;
    %     manually adding the mean value of data back among the modes
    out1.Phi = [Mnormalized(:), out1.Phi];
    out1.b = [bM; out1.b];
    out1.omega = [0; out1.omega];
    out1.lambda = [1; out1.lambda];
    %%
    Ur_data=out1.Phi;
    p = size(Ur_data,2);
    Nzeros=zeros(p,1);
elseif DataProjection == 3
    %AUS
    Ur_data=0;
    p=numModes;
    Nzeros=zeros(p,1);
end
end
