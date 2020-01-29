clear all;clc;
rng(1331); % set a random number seed for consistent simulations.
dimension = 100;  % dimension of the Lorenz96 system.
lorenzinit = rand(dimension,1);  % initial conditions.
outputtimes = linspace(0,10,dimension);  % output times for the ode45 calls.
[~,w] = ode45(@lorenz96,outputtimes,lorenzinit);
%% POD
tol=.0001;
lorenz96run = w';
% [r, Xr, Ur, Vr, Sr] = orderReduction(1-tol, lorenz96run); % Get POD
% save ('pod', 'Ur')
%% DMD
r_dmd=90;
[Phi,lambda] = orderReduction_DMD(r_dmd, lorenz96run);% Get DMD



% save ('DMD', 'Phi')