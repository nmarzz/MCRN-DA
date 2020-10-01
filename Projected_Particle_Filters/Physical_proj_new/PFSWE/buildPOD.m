function [Ur] = buildPOD(tolerance, model_output)
% model_output = model_output';
% [Ur] = orderReduction_POD(tolerance, model_output); % Get POD
[U,S,V] = svd(model_output,'econ');
Ur = U(:,1:tolerance);
% Sr = S(1:tolerance,1:tolerance);
% Vr = V(:,1:tolerance)'; % Truncate U,S,V using the rank r
% Xr = Ur*Sr*Vr; % Truncated matrix
end
