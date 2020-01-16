rng(1331); % set a random number seed for consistent simulations.
 
dimension = 400;  % dimension of the Lorenz96 system.
lorenzinit = rand(dimension,1);  % initial conditions.
outputtimes = linspace(0,10,400);  % output times for the ode45 calls.
[~,w] = ode45(@lorenz96,outputtimes,lorenzinit);
%% POD
tol=.001;
lorenz96run = w';
[r, Xr, Ur, Vr, Sr] = orderReduction(1-tol, lorenz96run); % Get POD

save ('pod', 'r', 'Xr','Ur','Vr','Sr')
 