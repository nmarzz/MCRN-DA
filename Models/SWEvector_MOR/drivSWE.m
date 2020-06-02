
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
forecast_length_days = 1;   % Total simulation length (days)

dt = dt_mins*60.0; % Timestep (s)
output_interval = output_interval_mins*60.0; % Time between outputs (s)
forecast_length = forecast_length_days*24.0*3600.0; % Forecast length (s)
nt = fix(forecast_length/dt)+1; % Number of timesteps
timesteps_between_outputs = fix(output_interval/dt);
noutput = ceil(nt/timesteps_between_outputs); % Number of output frames

% Specify the range of heights to plot in metres
plot_height_range = [9500 10500];

[XNEW,u,v,h,pars,parsanim]=init_swe;

nx=pars.nx;
ny=pars.ny;

dx=pars.dx;
dy=pars.dy;
x=parsanim.x;
y=parsanim.y;
H=parsanim.H;

%NUMBER OF LYAPUNOV EXPONENTS
p=1
%DIMENSION OF PROBLEM
DIM=3*nx*ny

%INITIALIZE
LE = zeros(p,1);
%q = eye(DIM,p);
q=orth(rand(DIM,p));
FEVAL = zeros(DIM,p);
NEWDIFF = zeros(DIM,p);
TIME = 0;
 

sqrteps = sqrt(eps(1));
%INTEGRATE SWE

T0=0;

i_save=1;

for i = 1:nt

   XOLD = XNEW;
%UPDATE TRAJECTORY
   XNEW = formod(TIME,XOLD,dt,pars);

% BEGIN: To animate
[u,v,h] = swe_vectortodata(XNEW, nx, ny);

x_save(:,i_save) = XNEW;

u_save(:,:,i_save) = u;
v_save(:,:,i_save) = v;
h_save(:,:,i_save) = h;
t_save(i_save) = (i-1).*dt;
i_save = i_save+1;
% END: To animate

% %% OVER ALL THE LYAPUNOV EXPONENTS WE WANT, SOME p<=DIM
%    for j=1:p
% %EVALUATE F(X+eps^{1/2}*Qj)
%       NEWIC = XOLD+sqrteps*q(:,j);
%       QTAU = formod(TIME,NEWIC,dt,pars);
%       NEWDIFF(:,j) = (QTAU - XNEW)/sqrteps;
%    end
% 
% %CALL mgs
%    [q,r] = mgs(NEWDIFF);
% 
% %FORM LES
%    for j=1:p
%        LE(j) = LE(j) + log(r(j,j));
%    end 
% 
% %UPDATE TIME 
%    TIME = TIME + dt

end

% LE = LE/(TIME-T0)

% save('SWE.mat','x_save','x','y','H','dt','t_save','plot_height_range')

% disp('Now run "animate_new" to animate the simulation');
