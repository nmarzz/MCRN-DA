function [Xa, W] = letkf(X, yobs, H, Rdiag, rho)
%LETKF - Local Ensemble Transform Kalman Filter.
%  This function implements the "cookbook" in the paper by
%  Brian R. Hunt, Eric J. Kostelich, and Istvan Szunyogh,
%  "Efficient data assimilation for spatiotemporal chaos:
%  A local ensemble transform Kalman filter," Physica D 230 (2007), 112-126.
%  Preprints may be found at <http://www.weatherchaos.umd.edu/publications.php>
%  or at <http://arxiv.org/abs/physics/0511236>.
%
%  Input arguments
%  ---------------
%  X : Initial background ensemble, represented as a matrix of size N x NSOL.
%      (NSOL is the number of ensemble solutions.)  In general, X is a portion
%      of the full model state; the localization is performed by the caller.
%  YOBS : A NOBS-vector of observations that the caller deems relevant to
%      updating the background X.  Required argument.  If YOBS is empty,
%      then LETKF returns the background, X.
%  H : The observation operator, represented as a NOBS x NSOL matrix.
%      More precisely, the Kth column of H represents the result of the
%      observation operator applied to X(:,k), i.e., to the Kth ensemble member.
%      Required argument.
%  RDIAG : NOBS-vector of errors in the observations Y.  In this implementation,
%      we assume that the observation errors are uncorrelated, i.e., that
%      the observation covariance matrix is diagonal.  See the paper for
%      appropriate adjustments for the non-diagonal case.  Required argument.
%  RHO : covariance inflation factor (scalar).  See the paper for further
%      discussion.  Optional argument; if no value is specified, then 
%      1 is assumed (i.e., no inflation is performed).
%
%  Output arguments
%  ----------------
%  Xa : The adjusted model state ("analysis") given the observations,
%       as a matrix of the same size as X.
%  W  : (optional)  The NSOL x NSOL weight matrix giving the required
%       linear combination of background ensemble perturbations needed
%       to form Xa.
%
%  Copyright 2009 by Eric J. Kostelich.
%  This program is free software: you may redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <http://www.gnu.org/licenses/>.

%  Begin LETKF.
    nobs = numel(yobs);  % = size(H, 1)
    nsol = size(X, 2);
    if(size(H, 1) ~= nobs || size(H, 2) ~= nsol)
       error('letkf: H must be NOBS x NSOL')
    end
    if(nargin < 5)  % then don't perform variance inflation
       rho = 1;
    end

%  If there are no observations, then return the background.
    if(nobs == 0)
       Xa = X;
       W = eye(nsol);  % identity matrix
       return
    end
 
%  Step 1.  Obtain the observation mean and perturbation matrix Yb.
    ybar = sum(H, 2) / nsol;  % column mean
    Yb = H;
    for k = 1:nsol
       Yb(:,k) = Yb(:,k) - ybar;
    end

%  Step 2.  Obtain the background mean and perturbation matrix Xb.
    xbar = sum(X, 2) / nsol;  % column mean
    Xb = X;
    for k = 1:nsol
       Xb(:,k) = Xb(:,k) - xbar;
    end

%  Step 3.  Compute Yb and assume that R is diagonal to form R^(-1)*Yb.
    RinvYb = Yb;
    for k = 1:nsol
       RinvYb(:,k) = RinvYb(:,k) ./ Rdiag; 
    end

%  Step 4.  Deferred to below.
%  Step 5A.  Compute the NSOL x NSOL matrix Q = [(k-1)I/rho + C*Yb].
%  Since Q is symmetric, its eigendecomposition is Q = V D V', which
%  can be used efficiently to form its inverse and symmetric square root.
    Q = Yb' * RinvYb;
    inflation = (nsol - 1) / rho;  % variance inflation, if any
    for k = 1:nsol
       Q(k,k) = Q(k,k) + inflation;
    end
    [V, D] = eig(Q);  %  V is the unitary matrix of eigenvectors of Q

%  Step 5B.  Form Pa~ = inv(Q).  It's possible that small negative
%  eigenvalues can be created due to roundoff, so set them to zero.
    lambda = max(diag(D), 0);
    Z = zeros(nsol, nsol);
    for k = 1:nsol
       if(lambda(k) > 0)
          Z(:,k) = V(:,k) / lambda(k);
       end
    end
    Pa = Z * V';

%  Step 6.  Compute W = sqrt((k-1)*Pa~).
    for k = 1:nsol
       if(lambda(k) > 0)
          Z(:,k) = V(:,k) / sqrt(lambda(k));
       else
          Z(:,k) = 0;  % just in case
       end
    end
    W = (Z * V') * sqrt(nsol-1);

%  Steps 7 and 4.  Form wbar and complete the computation of W.
%  Here we assume a diagonal observation covariance matrix.
%  See the paper for appropriate adjustements when the observation
%  errors are not independent.
    wbar = Pa * (Yb' * ((yobs - ybar) ./ Rdiag));  % elementwise division
    for k = 1:nsol
       W(:,k) = W(:,k) + wbar;
    end
   
%  Step 8.  Form the analysis ensemble.
    Xa = Xb * W;
    for k = 1:nsol
       Xa(:,k) = Xa(:,k) + xbar;
    end
    return
end  % letkf
