function out = normlogw( x )
% NORMLOGW Normalize a vector of logarithmic variables so that 
% sum(exp(x)) = 1
%
% I'm trying to compute
% 
% w = exp(x);
% w = w/sum(w);
% outx = log(w);
%
% without NaN-overflow.
%
% This "trick" normalizes all weights by the largest of them
% (x - max(x)) is normalization in log-domain
% and then performs the computation so that the max value that is
% exponentiated is always just 0 ( exp( x-max(x) ) <= 1 .

X = max(x);
xnorm = x - X;
out = xnorm - log(sum(exp(xnorm)));

