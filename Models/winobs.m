function [y, R, H] = winobs(yobs, Rdiag, choice, background, window)
%WINOBS - Determine the subset of the observations that lies in the current
%  local window.  Internal function.

   [common, ic, iw] = intersect(choice, window);  % MATLAB builtin
   y = yobs(ic);
   common = choice(ic);  % so they match the order in Y
   H = background(common,:);
   R = Rdiag(ic);
   return
end  % winobs
