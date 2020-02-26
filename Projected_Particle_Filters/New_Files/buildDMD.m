function [Phi]=buildDMD(numModes,model_output,dt)
% DMD
model_output = model_output';
[Phi] = orderReduction_DMD(numModes, model_output,dt);% Get DMD
end
