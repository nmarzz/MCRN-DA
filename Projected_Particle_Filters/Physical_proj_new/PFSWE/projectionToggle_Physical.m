function [V] = projectionToggle_Physical(PhysicalModelProjection,Ur,p)
% Type of projection for the Physical Model
% V: orthonormal columns, physical model 
% p: rank of projection
if PhysicalModelProjection == 0
%     V = 1;
      V = sparse(eye(p));
elseif PhysicalModelProjection == 1
    V = getpod(Ur,p);
elseif PhysicalModelProjection == 2
    V = getDMD( Ur,p );
end

assert(isreal(V),"Projection matrix has to be real-valued");
