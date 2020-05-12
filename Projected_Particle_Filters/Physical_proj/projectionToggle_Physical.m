function [V] = projectionToggle_Physical(PhysicalModelProjection,Model_Dimension,Ur,p)
% Type of projection for the Physical Model
% V: orthonormal columns, physical model 
% p: rank of projection
if PhysicalModelProjection == 0
    V = eye(Model_Dimension);   
elseif PhysicalModelProjection == 1
    V = getpod(Ur,p);
elseif PhysicalModelProjection == 2
    V = getDMD( Ur,p );
end

