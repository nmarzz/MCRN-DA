function [Phi]=buildDMD(numModes,model_output)
% DMD
model_output = model_output';
[Phi] = orderReduction_DMD(numModes, model_output);% Get DMD
end
