function [analysis,xbar,background]=letkfstep(F,time,h,irad,background,yobs,choice,Rdiag,rho)
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

       % Generate analysis
       [analysis,xbar] = create_analysis(irad, background, yobs, choice, Rdiag, rho);

       % Generate the next forecast
 %      background = rkfixed(F, time, analysis, h);
 sz = size(background);
 for ii = 1:sz(2)
     background(:,ii) = formod(time, analysis(:,ii), h, F);
 end
