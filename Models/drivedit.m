% Driver file to set up the system, assign its dimensions, and call the run function.

[truth,u,v,h,pars,parsanim]=init_swe; % Initialize SWE
F = pars;       % I'll pass  a bunch of parameters, held in pars, on to run.
% F=@lorenz96;  % This was for Lorenz96.
dim=38100;  % number of dimensions.
irad=5;  % I assume irad is iterates.
runnew(F,truth,dim,20,200,irad,0.5,1,1)

