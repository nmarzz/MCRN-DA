function [Ur_physical,p,Nzeros] = Projection_physical_type(PhysicalProjection,numModes,tolerance,Model_Dimension,Built_Model,dt)

%Physical Projection
if PhysicalProjection == 0
    Ur_physical=0;
    p=Model_Dimension;
    Nzeros=zeros(Model_Dimension,1);
elseif PhysicalProjection == 1
    %POD
    Ur_physical= buildPOD(tolerance,Built_Model);
    p = size(Ur_physical,2);
    Nzeros=zeros(p,1);
elseif PhysicalProjection == 2
    %DMD
    DataMatrix=Built_Model;%snapshot matrix
    r=numModes;
    VALUE=true;
    M = mean(DataMatrix, 2); % 2 = mean across columns
    stepVALUE=-1;
    out1 = dmd( DataMatrix, dt, r, 'removemean', VALUE,'step',stepVALUE,'sortbyb', VALUE);
% %     out = dmd( DataMatrix, dt, r);
    bM = norm(M);
    Mnormalized = M/bM;
    % manually adding the mean value of data back among the modes
    out1.Phi = [Mnormalized(:), out1.Phi];
    out1.b = [bM; out1.b];
    out1.omega = [0; out1.omega];
    out1.lambda = [1; out1.lambda];
    
    Ur_physical=out1.Phi;
    p = size(Ur_physical,2);
    Nzeros=zeros(p,1);
    
elseif PhysicalProjection == 3
    %%AUS
    Ur_physical=0;
    p=numModes;
    Nzeros=zeros(p,1);
end
end

