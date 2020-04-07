function xp = lorenz96(t, x)
% Adding a test comment
%LORENZ96 - The vector field for the Lorenz '96 system.
%  Input arguments
%  ---------------
%  X : the current state, generally an ensemble of solutions stored as an
%      N x NSOL matrix.  (N is the system size of the Lorenz '96 system
%      and NSOL is the ensemble size.)
%  T : the current time.  Although the system is autonomous, we include the
%      time for compatibility with general ODE solvers.
%
%  Ouput argument
%  --------------
%  XP - an N x NSOL matrix whose Kth column is the time rate of change
%       associated with the Kth ensemble solution.
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

   n = size(x, 1);  % typically 40 but could be something else
   F = 8;   % forcing, fixed here to Lorenz's value
   jp1 = circshift([1:n]', -1);
   jm1 = circshift([1:n]', 1);
   jm2 = circshift([1:n]', 2);
   xp = (x(jp1,:) - x(jm2,:)) .* x(jm1,:) - x + F;
   return
end  % lorenz96
