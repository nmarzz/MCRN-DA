% Driver file to set up the system, assign its dimensions, and call the run function.


% F=@lorenz96;  % Needs to setup something about SWE to indicate we're using that.
dim=40;  % number of dimensions sounds find.
irad=5;  % I assume irad is iterates.
run(F,dim,20,200,irad,0.5,1,1)

