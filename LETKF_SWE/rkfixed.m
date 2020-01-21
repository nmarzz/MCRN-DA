function y = rkfixed(f, t, x, h)
%RKFIXED - Perform one step of a fixed-step, fourth-order Runge-Kutta method.
%  Input arguments
%  --------------
%  F : the vector field (as a function handle)
%  T : current time
%  X : current solution.  An ensemble of solutions is supported: each column
%      of X is assumed to be one current solution, and the RK method is
%      applied columnwise.
%  H : the time step.
%  Output argument
%  ---------------
%  Y : The solution X(T) advanced to X(T+H) under the vector field F.
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

%  Begin RKFIXED.
    k1 = f(t, x);
    k2 = f(t + (h/2), x + (h/2)*k1);
    k3 = f(t + (h/2), x + (h/2)*k2);
    k4 = f(t + h, x + h*k3);
    y = x + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
    return
end  % rkfixed
