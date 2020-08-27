function Y = lorenz96(x0, T, F, eta, P)
% LORENZ96 Given initial condition of N elements, forcing magnitude F, and
% time vector, system simulates the Lorenz'96 model with periodic forcing (not to be
% confused with the more-famous Lorenz'63).
%
% Lorenz'96 would be the model where F is constant.
% Here we replace constant F with F*(1+sin(2*pi*t/P)/2)
% Therefore the forcins oscillates between (F/2 and 3F/2)
%
% Good choices for length of x0 are N = 40, 55, 60
% Good choices for F are 3. 3.5, 4, 8
% P - likely doesn't do much, but can be made to turn chaos on and off.
%
% After simulation, normal zero-mean noise of magnitude eta is added.
%
% Returns X - numel(x0) x numel(T) matrix of snapshot evolution.


N = numel(x0);
omega = 2*pi/P;

% simulate deterministic system
[~,X] = ode45( @(t,x)lorenz96_rhs(t,x,F,omega), T, x0(:), ...
    odeset('Vectorized','on','AbsTol',1e-10, 'RelTol',1e-6));

X = transpose(X);
Y = X + randn(size(X))*eta; % add noise

% if no outputs were requested, visualize.
if nargout == 0
    tiledlayout('flow');
    xax = transpose(1:N);
    nexttile; pcolor(T,xax, Y); shading flat; colorbar;
    nexttile; waterfall(xax,T, Y'); shading flat; colorbar;
    
    clear X;
end

end

function dx = lorenz96_rhs(t,x,F,omega)

dx = ( circshift(x,1,1) - circshift(x,-2,1) ).*circshift(x,-1,1) - x + F*(2+sin(omega*t))/2;

end
