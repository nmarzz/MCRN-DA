function [ensemble] = createbackground(truth, pert, nsol)
%INITIALIZE - initialize an ensemble givne a "truth".
%  The ensemble is similar to the truth, except that
%  additional Gaussian noise is added.  
%
%  Input arguments
%  ---------------
%  PERT : the perturbation size
%  NSOL : the ensemble size.
%  TRUTH : the "true" state for comparison with an analysis/forecast system.
%  Output arguments
%  ----------------
%  ENSEMBLE : the background ensemble 

 ensemble = repmat(truth,[1,nsol]);
 ensemble = ensemble + pert*randn(size(ensemble));
     
 end 
