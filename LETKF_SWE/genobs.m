function [yobs, choice, Rdiag] = genobs(x, obsfrac, variance)
%GENOBS - generate synthetic observations of the Lorenz96 system.
%  We pick random locations along the "truth", X, that are observed.
%  The observations are simply the truth plus a Gaussian random value
%  of mean 0 with the given VARIANCE.
%  CHOICE holds the chosen locations (and allows us to implement the
%  observation operator).
%  RDIAG is simply the diagonal of the observation covariance matrix
%    (the observation errors are independent).
%  Internal function.
%
   n = size(x, 1);   % size of Lorenz system
   p = min(n, max(0, round(n*obsfrac)));
   c = randperm(n);
   choice = c(1:p);
   yobs = x(choice) + randn(p, 1) * sqrt(variance);
   Rdiag = repmat(variance, size(yobs));
   return
end  % genobs
