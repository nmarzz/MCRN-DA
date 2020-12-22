function [out1,filename] = simulationDMD(numModes,Built_Model,dt)
% this is something you'd write to get out1.
% >> [out1] = simulationPOD(tolerance,model_output);
% However, every time
% you run this script, you have to wait 5 seconds.
% this is a better way - let's pick a file name based on parameters:
filename = sprintf("DMD_r%02d.mat", numModes);

if exist(filename,'file') % if file exists, load out1
    load(filename)
    disp("Loaded from "+filename);
else
    DataMatrix=Built_Model;%snapshot matrix
    Nsteps = size(DataMatrix, 2); % number of columns
    r=numModes;
    VALUE=true;
    M = mean(DataMatrix, 2); % 2 = mean across columns
    stepVALUE = round( linspace(1, Nsteps, 5) ); % steps at which b coefficients will be computed
    %stepVALUE=[1 100 200 300 500];
%     out1 = dmd( DataMatrix, dt, r, 'removemean', VALUE, 'sortbyb', VALUE);
    out1 = dmd( DataMatrix, dt, r-1, 'removemean', VALUE,'step',stepVALUE,'sortbyb', VALUE);
    bM = norm(M);
    Mnormalized = M/bM;
    % manually adding the mean value of data back among the modes
    out1.Phi = [Mnormalized(:), out1.Phi];
    out1.b = [bM; out1.b];
    out1.omega = [0; out1.omega];
    out1.lambda = [1; out1.lambda];

    disp('Yep, simulating-DMD hard!');
    save(filename,'out1') % store into filename our output variables
    disp("Saved to "+filename);
end

end
% close all; scatter(real(out1.lambda), imag(out1.lambda),1+50*abs(out1.b(:,2))/max(abs(out1.b(:,2))),'filled')