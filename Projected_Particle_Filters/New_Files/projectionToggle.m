function [q_data,q_physical] = projectionToggle(PhysicalModelProjection,DataModelProjection,Model_Dimension,Ur,p,Phi)
% Type of projection for the Data Models
if DataModelProjection == 0
    q_data = eye(Model_Dimension); %identity matrix
elseif DataModelProjection == 1
    q_data = getpod(Ur,p);
elseif DataModelProjection == 2
    q_data = getDMD( Phi,p );
elseif DataModelProjection == 3
    %then q is the AUS basis
end

% Type of projection for the Physical Model
if PhysicalModelProjection == 0
    q_physical = eye(Model_Dimension);
elseif PhysicalModelProjection == 1
    q_data = getpod(Ur,p);
elseif PhysicalModelProjection == 2
    q_data = getDMD( Phi,p );
end
