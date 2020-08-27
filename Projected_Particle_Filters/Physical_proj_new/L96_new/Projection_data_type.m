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
    %numModes=300;  % number of DMD_modes you want to use
    Ur_data=buildDMD(numModes,Built_Model,dt);
    p = size(Ur_data,2);
    Nzeros=zeros(p,1);
elseif DataProjection == 3
    %AUS
    Ur_data=0;
    p=numModes;
    Nzeros=zeros(p,1);
end
end
