function [q_data,q_physical] = projectionToggle(PhysicalModelProjection,DataModelProjection,Model_Dimension)
% Type of projection for the Data Models
if DataModelProjection == 0
    q_data = eye(Model_Dimension); %identity matrix
elseif DataModelProjection == 1
    %then q is POD basis
elseif DataModelProjection == 2
    %then q is the DMD basis
elseif DataModelProjection == 3
    %then q is the AUS basis
end

% Type of projection for the Physical Model
if PhysicalModelProjection == 0
    q_physical = eye(Model_Dimension); 
elseif PhysicalModelProjection == 1
    %then q is POD basis
elseif PhysicalModelProjection == 2
    %then q is the DMD basis
end

