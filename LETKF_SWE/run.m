function [RMSE_save] = run(xnew,pars,parsanim,dim,nsol, nsteps, irad, obsfrac, obserr, moderr, rho)
%RUN96 - Run an analysis/forecast cycle for the Lorenz96 model.
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

%%% UPDATE: 
    [truth, background] = initialize(xnew,pars,dim,nsol);
    figure(1);
    LEs = zeros(nsol,1);
    h=0.05;
    xvals = linspace(1,dim,dim);

    for k = 1:nsteps
       time=k*h;

       % Generate synthetic observations from the truth and assimilate them
       % NEED: fix localization for a 2D spatial domain.
       [yobs, choice, Rdiag] = genobs(truth, obsfrac, obserr);

       % LETKF Step: NEED to update to add in model error.
       [analysis,xbar,background]=letkfstep(pars,time,h,irad,background,yobs,choice,Rdiag,moderr,rho);

%%%%%% REPLACED with call to formod
%       truth = rkfixed(F, time, truth, h);
        truth = formod(time,truth,h,pars);
        RMSE = norm(truth-xbar,2)/sqrt(dim);
        RMSE_save(k) = RMSE;

       [q,r]=mgs(analysis); 
       LEs = LEs + log(diag(r));

       animate(time,truth-real(xbar),pars,parsanim);
       %animate(time,real(xbar),pars,parsanim);
       makeplot(nsol,xvals,truth,xbar,analysis,choice,yobs,k,h,RMSE)
       pause(0.1)

    LocalLEs = LEs/(k*h)
    end

   LEs = LEs/(nsteps*h)
   return
end
