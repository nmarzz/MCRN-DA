
function [q_data] = projectionToggle_data(DataModelProjection,Model_Dimension,Ur,p)
% Type of projection for the Data Models
% U: orthonormal columns,data model
% p: rank of projection
if DataModelProjection == 0
    U = eye(Model_Dimension); %identity matrix
    q_data = eye(Model_Dimension); %identity matrix
elseif DataModelProjection == 1
    U = getpod(Ur,p);
elseif DataModelProjection == 2
    U = getDMD( Ur,p );
    q_data = getDMD( Ur,p );
    % elseif DataModelProjection == 3
    %    [U,~] = getausproj(Model_Dimension,p,Fmod,t,est,h,q,LE);
end


