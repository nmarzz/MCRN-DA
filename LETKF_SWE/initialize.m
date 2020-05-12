function [truth, ensemble] = initialize(xnew,pars, n, nsol)
%INITIALIZE - initialize an ensemble and a corresponding "truth" for the
%  Lorenz '96 system.  The truth is generated from a small pertubation of
%  the unstable equilibrium points; the ensemble is similar, except that
%  additional Gaussian noise is added.  The integration is for 500 time
%  steps of size 0.05 (as in Lorenz's paper), corresponding to 125 days
%  of "weather".
%
%  Input arguments
%  ---------------
%  N : the size of the Lorenz '96 system (typically 40, but can be another
%      value.
%  NSOL : the ensemble size.
%  Output arguments
%  ----------------
%  TRUTH : the "true" state for comparison with an analysis/forecast system.
%  ENSEMBLE : the background ensemble (essentially a climatology for the
%    Lorenz '96 system).

%Lorenz 96
%    forcing = 8;  % Lorenz's value
%    pert = forcing / 10;  % something small

     pert = 1;

%  Start with a small pertubation at the middle point from the equilibrium
%    ensemble = repmat(forcing, n, nsol);
%    truth = ensemble(:,1);
%    truth(n/2) = truth(n/2) + pert;  % middle point
%    ensemble(n/2,:) = ensemble(n/2,:) + pert * randn(1, nsol);
   
%Lorenz 63
%    ensemble = repmat(forcing, n, nsol);
%    truth = ensemble(:,1);
%    ensemble(2,:) = ensemble(2,:) + pert * randn(1, nsol);

%Shallow Water Eqnts.  SWE.  This doesn't currently work.
     truth=xnew;
     ensemble = repmat(truth,[1,nsol]);
     ensemble = ensemble + pert*randn(size(ensemble));

%  Spin up
%    h = 0.05;  % Lorenz's value
%    for k = 1:500
%       truth = rkfixed(@lorenz96, k * 0.05, truth, h);
%       ensemble = rkfixed(@lorenz96, k * 0.05, ensemble, h);
%    end
    return
end  % initialize
