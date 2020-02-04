function buildDMD(numModes,dimension,model)

% Dimension is in input to function

rng(1331); % set a random number seed for consistent simulations.
lorenzinit = rand(dimension,1);  % initial conditions.
outputtimes = linspace(0,10,dimension);  % output times for the ode45 calls.
[~,w] = ode45(model,outputtimes,lorenzinit); % Can change the model

lorenz96run = w';

% DMD
[Phi] = orderReduction_DMD(numModes, lorenz96run);% Get DMD
save ('DMD', 'Phi')
end

