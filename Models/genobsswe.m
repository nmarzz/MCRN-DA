% Comment from runSWE.m:
%  Generate synthetic observations from the truth and assimilate them
% [yobs, choice, Rdiag] = genobsSWE(truth, obsfrac, obserr);

function [yobs, choice, Rdiag] = genobsSWE(x, obsfrac, variance)
%GENOBS - generate synthetic observations of the shallow water equations
%  We pick random locations along the "truth", X, that are observed.
%  The observations are simply the truth plus a Gaussian random value
%  of mean 0 with the given VARIANCE.
%  CHOICE holds the chosen locations (and allows us to implement the
%  observation operator).
%  RDIAG is simply the diagonal of the observation covariance matrix
%    (the observation errors are independent).
%  Internal function.
%
   n = size(x, 1);   % size of system
   p = min(n, max(0, round(n*obsfrac)));
   c = randperm(n);
   choice = c(1:p);
   yobs = x(choice) + randn(p, 1) * sqrt(variance);
   Rdiag = repmat(variance, size(yobs));
   return
end  % genobs
