
% Copyright (c) 2014 by Robin Hogan
%
% Copying and distribution of this file, with or without modification,
% are permitted in any medium without royalty provided the copyright
% notice and this notice are preserved.  This file is offered as-is,
% without any warranty.
%
% This model integrates the shallow water equations in conservative form
% in a channel using the Lax-Wendroff method.  It can be used to
% illustrate a number of meteorological phenomena.

dt_mins              = 1;   % Timestep (minutes)
output_interval_mins = 60;  % Time between outputs (minutes)
forecast_length_days = 4;   % Total simulation length (days)

dt = dt_mins*60.0; % Timestep (s)
output_interval = output_interval_mins*60.0; % Time between outputs (s)
forecast_length = forecast_length_days*24.0*3600.0; % Forecast length (s)
nt = fix(forecast_length/dt)+1; % Number of timesteps
timesteps_between_outputs = fix(output_interval/dt);
noutput = ceil(nt/timesteps_between_outputs); % Number of output frames

% Specify the range of heights to plot in metres
plot_height_range = [9500 10500];

%Initialize: IC, parameters, and parameters for visualizing
[XNEW,u,v,h,pars,parsanim]=init_swe;
%loop to solve the DE then POD

%Spatial discretization
nx=pars.nx;
ny=pars.ny;
dx=pars.dx;
dy=pars.dy;

%Parameters for animation
x=parsanim.x;
y=parsanim.y;
H=parsanim.H;

%DIMENSION OF PROBLEM
DIM=3*nx*ny

irad=4; %Currently need 2*irad+1 <=min(nx,ny)
%Run LETKF
%function run(xnew,pars,parsanim,dim,nsol, nsteps, irad, obsfrac, obserr, moderr, rho)
%rho<=1
RMSE_save = run(XNEW,pars,parsanim,DIM,40,50,irad,0.1,1.0,1.0,1)
%Requires a routine to time step SWE of the form:
%function xnew= formod(t,x,dt,pars)

figure(2);
plot(RMSE_save);

%Move animation stuff into run so that we can visualize
%or add here for post simulation visualization.
