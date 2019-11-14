function run(F,dim,nsol, nsteps, irad, obsfrac, obserr, rho)
%  RUN - Run an analysis/forecast cycle for the model.
%  NSOL : the number of forecast ensembles to create and analyze.
%  NSTEPS : the number of time steps that the cycle should run.
%     (One step corresponds to about 6 hours of "weather".)
%  OBSFRAC : the fraction of points at which observations are taken.
%  OBSERR : the observation error (variance) in the generated observations.
%  RHO : the covariance inflation factor applied to the analysis.
%  This function generates two graphs: the first shows the truth,
%  analysis mean, and generated observations, and the second shows
%  the individual ensemble analyses.  Each graph is updated at each step
%  of the cycle.
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

    [truth, background] = initswe(F,dim,nsol);
    %h = 0.05;
    figure(1);
    %pause
    LEs = zeros(nsol,1);
    h=0.05;
    xvals = linspace(1,dim,dim);

    for k = 1:nsteps
       time=k*h;
       %  Generate synthetic observations from the truth and assimilate them
       [yobs, choice, Rdiag] = genobsswe(truth, obsfrac, obserr);

       [analysis,xbar,background]=letkfstep(F,time,h,irad,background,yobs,choice,Rdiag,rho);

       %[analysis,xbar] = create_analysis(irad, background, yobs, choice, Rdiag, rho);
       %% Generate the next forecast
       %h=0.05;
       %background = rkfixed(F, k * 0.05, analysis, h);

       truth = rkfixed(F, time, truth, h);

       [q,r]=mgs(analysis); 
       LEs = LEs + log(diag(r));

       %  Plot the truth, analysis mean, and observations on graph #1
       subplot(1,2,1);
       plot(truth, 'k-');
       hold on
          %axis([1 40 -10 15]);
          plot(xvals,xbar, 'b--');
          %size(xbar)
          %pause
          plot(choice, yobs, 'r*');
          title(['Time = ',num2str(k*h)])
       hold off

       %  Plot the ensemble and observations on graph #2
       subplot(1,2,2); 
       plot(xvals,analysis(:,1), 'b-');
       hold on
          %axis([1 40 -10 15]);
          %size(analysis)
          %pause
          for j = 2:nsol
             plot(xvals,analysis(:,j), 'b-');
          end
          plot(choice, yobs, 'r*');
          title(['Time = ',num2str(k*h)])
       hold off
        pause(0.1)
%       if(k == 1)
%          pause
%       else
%          pause(0.1)
%       end

%       % Generate the next forecast
%       h=0.05;
%       truth = rkfixed(F, k * 0.05, truth, h);
%       background = rkfixed(F, k * 0.05, analysis, h);
    end

   LEs = LEs/(nsteps*h)
    return
end
