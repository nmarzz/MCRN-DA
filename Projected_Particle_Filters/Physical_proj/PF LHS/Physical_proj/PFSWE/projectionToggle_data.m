function [U] = projectionToggle_data(DataModelProjection,Ur,p)
% Type of projection for the Data Models
% U: orthonormal columns,data model
% p: rank of projection
if DataModelProjection == 0
    U = 1;
elseif DataModelProjection == 1
    U = getpod(Ur,p);
elseif DataModelProjection == 2
    U = getDMD( Ur,p );
% elseif DataModelProjection == 3
%    [U,~] = getausproj(Model_Dimension,p,Fmod,t,est,h,q,LE);
end


