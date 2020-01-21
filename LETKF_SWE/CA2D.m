function [Xa,xbar] = CA2D(irad, background, yobs, choice, Rdiag, rho,nx,ny)
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

%To do: iradx and irady

   wlen = 2*irad + 1;
   n = size(background, 1);
   nsol = size(background, 2);
   Xa = zeros(size(background));

   for i=1:3
   for j = 1:nx
   for k = 1:ny

      [windowj,pointj]=winpt(j,irad,wlen,nx);
      [windowk,pointk]=winpt(k,irad,wlen,ny);

      point=(i-1)*nx*ny+mattovec(pointj,pointk,nx,ny);

      window(1)=(i-1)*nx*ny + mattovec(j,k,nx,ny);
      winjk=window(1);    
      ii=2;
      for l=1:wlen
          val=(i-1)*nx*ny + mattovec(j,windowk(l),nx,ny); 
          if val ~= winjk
             window(ii)=val;
             ii=ii+1;
          end
          %
          val=(i-1)*nx*ny + mattovec(windowj(l),k,nx,ny); 
          if val ~= winjk
             window(ii)=val;
             ii=ii+1;
          end
      end

      [y, R, H] = winobs(yobs, Rdiag, choice, background, window);
      analysis = letkf(background(window,:), y, H, R, rho);
      kk=(i-1)*nx*ny+mattovec(j,k,nx,ny);
      %point=wlen; %Always use midpoint
      point=1; %Since window(1) is where we store the current point.
      Xa(kk,:) = analysis(point,:);  % center point of the window.
   end
   end
   end
   xbar = sum(Xa, 2) / nsol;
   return
end  % create_analysis
