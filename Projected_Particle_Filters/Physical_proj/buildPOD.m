function [Ur] = buildPOD(tolerance, model_output)
model_output = model_output';
[Ur] = orderReduction_POD(tolerance, model_output); % Get POD
end
