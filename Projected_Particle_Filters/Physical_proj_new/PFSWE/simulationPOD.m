function [out1,filename] = simulationPOD(tolerance,model_output)
% this is something you'd write to get out1.
% >> [out1] = simulationPOD(tolerance,model_output);
% However, every time
% you run this script, you have to wait 5 seconds.
% this is a better way - let's pick a file name based on parameters:
filename = sprintf("POD_r%02d.mat", tolerance);

if exist(filename,'file') % if file exists, load out1
    load(filename)
    disp("Loaded from "+filename);
else
    [out1] =buildPOD(tolerance, model_output); % Get POD
    disp('Yep, simulating-POD hard!');
    save(filename,'out1') % store into filename our output variables
    disp("Saved to "+filename);
end

end