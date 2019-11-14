function [Xa,xbar] = create_analysis(irad, background, yobs, choice, Rdiag, rho)
%CREATE_ANALYSIS - compute an analysis for the Lorenz96 model.
%  Input arguments
%  ---------------
%  IRAD : the number of points ("radius") on each size of the point to be
%    analyzed.  That is, the Kth spatial point is analyzed using points
%    K-IRAD:K+IRAD; the code handles the wraparound (point N+1 <=> point 1).
%  BACKGROUND : the N x NSOL background ensemble solutions.  There are NSOL
%    members of the ensemble, and each ensemble represents the N states of
%    the Lorenz96 model.
%  RDIAG : We assume uncorrelated observation errors, so RDIAG(k) is the
%    error (variance) in the Kth observation.
%  RHO : Variance inflation factor (use 1 for no inflation).
%  Output argument
%  ---------------
%  XA : The updated solution (analysis), computed pointwise for each
%       background ensemble solution.
%  xbar : analysis mean

   wlen = 2*irad + 1;
   n = size(background, 1);
   nsol = size(background, 2);
   Xa = zeros(size(background));
   for k = 1:n
%This depends on having periodic BCs
%      window = circshift([1:n]', irad-k+1);
%      point = irad+1;
%This works with non-periodic BCs
      if k-irad >= 1 & k+irad <=n
         window = circshift([1:n]', irad-k+1); %will use center point
         point = irad+1;
      elseif k+irad <=n
         window = [1:n]'; %will use "kth" point
         point = k;
      else 
         window = [n-2*irad:n]'; %will use kth point
         point=wlen+k-n;
      end
%End: non-periodic BCs
      window = window(1:wlen);  % define the local region
      [y, R, H] = winobs(yobs, Rdiag, choice, background, window);
      analysis = letkf(background(window,:), y, H, R, rho);
      Xa(k,:) = analysis(point,:);  % center point of the window.
   end
   xbar = sum(Xa, 2) / nsol;
   return
end  % create_analysis
