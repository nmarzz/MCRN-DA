function [analysis,xbar,background]=letkfstep(pars,time,h,irad,background,yobs,choice,Rdiag,moderr,rho)
% Input:
%  F - RHS of state space model
%  time - current time
%  h - timestep
%  IRAD : the number of points ("radius") on each size of the point to be
%    analyzed.  That is, the Kth spatial point is analyzed using points
%    K-IRAD:K+IRAD; the code handles the wraparound (point N+1 <=> point 1).
%  BACKGROUND : the N x NSOL background ensemble solutions.  There are NSOL
%    members of the ensemble, and each ensemble represents the N states of
%    the Lorenz96 model.
%  yobs - current observations
%  CHOICE holds the chosen locations (and allows us to implement the
%  observation operator).
%  RDIAG : We assume uncorrelated observation errors, so RDIAG(k) is the
%    error (variance) in the Kth observation.
%  RHO : Variance inflation factor (use 1 for no inflation).
%
% Output:
% analysis - the ensemble analysis
% xbar - the mean of the ensemble analysis
% background - the updated background
       nx=pars.nx;
       ny=pars.ny;
       [dim,nsol]=size(background);

       % Generate analysis
       %[analysis,xbar] = create_analysis(irad, background, yobs, choice, Rdiag, rho);
       %[analysis,xbar] = CA2D(irad, background, yobs, choice, Rdiag, rho, nx, ny);
       [analysis,xbar] = CA2D2(irad, background, yobs, choice, Rdiag, rho, nx, ny);

       % Generate the next forecast
%       background = rkfixed(F, time, analysis, h);
%      REPLACE with formod
       for j=1:nsol
       background(:,j) = formod(time,analysis(:,j),h,pars); % + mvnrnd(zeros(1,dim),moderr*eye(dim))'; % + MODEL ERROR;
       end
       %ADD MODEL ERROR
       %Slower:
       %background = background + mvnrnd(zeros(1,dim),moderr*eye(dim),nsol)';
       %Faster (only works with diagonal model error covariance matrix):
       background = background + mvnrnd(zeros(1,nsol),moderr*eye(nsol),dim);
%%%%%% REPLACE with call to formod
