function [w] = buildModel(dimension,model,Numsteps,T)
% Dimension is in input to function
rng(1331); % set a random number seed for consistent simulations.
% lorenzinit=cos( 2*pi*(1:dimension)/dimension);
lorenzinit = rand(dimension,1);  % initial conditions.
outputtimes = linspace(0,T,Numsteps);% output times for the ode45 calls.
[~,w] = ode45(model,outputtimes,lorenzinit); % Can change the model
end

