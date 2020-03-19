function [w] = buildModel(dimension,model,dt)
% Dimension is in input to function
rng(1331); % set a random number seed for consistent simulations.
lorenzinit = rand(dimension,1);  % initial conditions.
outputtimes = linspace(0,dt,dimension);% output times for the ode45 calls.
[~,w] = ode45(model,outputtimes,lorenzinit); % Can change the model
end

