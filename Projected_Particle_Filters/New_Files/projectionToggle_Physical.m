function [q_physical] = projectionToggle_Physical(PhysicalModelProjection,Model_Dimension,Ur,p)
% Type of projection for the Physical Model
if PhysicalModelProjection == 0
    q_physical = eye(Model_Dimension);
elseif PhysicalModelProjection == 1
    q_physical = getpod(Ur,p);
elseif PhysicalModelProjection == 2
    q_physical = getDMD( Ur,p );
end

