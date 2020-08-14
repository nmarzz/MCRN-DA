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
    %numModes=300;  % number of DMD_modes you want to use
    Ur_physical=buildDMD(numModes,Built_Model,dt);
    p = size(Ur_physical,2);
    Nzeros=zeros(p,1);
elseif PhysicalProjection == 3
    %AUS
    Ur_physical=0;
    p=numModes;
    Nzeros=zeros(p,1);
end
end

