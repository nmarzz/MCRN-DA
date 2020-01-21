function xnew= formod(t,x,dt,pars)

F=pars.F;
g=pars.g;
dx=pars.dx;
dy=pars.dy;
H=pars.H;
nx=pars.nx;
ny=pars.ny;

[u,v,h]=xtouvh(x,nx,ny);
%%Break up x into u,v,h
%u1=x(1:nx*ny);
%v1=x(nx*ny+1:2*nx*ny);
%h1=x(2*nx*ny+1:3*nx*ny);
%u=reshape(u1,nx,ny);
%v=reshape(v1,nx,ny);
%h=reshape(h1,nx,ny);

  % Compute the accelerations
  u_accel = F(2:end-1,2:end-1).*v(2:end-1,2:end-1) ...
              - (g/(2*dx)).*(H(3:end,2:end-1)-H(1:end-2,2:end-1));
  v_accel = -F(2:end-1,2:end-1).*u(2:end-1,2:end-1) ...
              - (g/(2*dy)).*(H(2:end-1,3:end)-H(2:end-1,1:end-2));

  % Call the Lax-Wendroff scheme to move forward one timestep
  [unew, vnew, h_new] = lax_wendroff(dx, dy, dt, g, u, v, h, ...
                                     u_accel, v_accel);

  % Update the wind and height fields, taking care to enforce
  % boundary conditions
  u = unew([end 1:end 1],[1 1:end end]);
  v = vnew([end 1:end 1],[1 1:end end]);
  v(:,[1 end]) = 0;
  h(:,2:end-1) = h_new([end 1:end 1],:);

xnew = uvhtox(u,v,h,nx,ny);
%%Put into xnew
%u1=reshape(u,nx*ny,1);
%v1=reshape(v,nx*ny,1);
%h1=reshape(h,nx*ny,1);
%xnew=cat(1,u1,v1,h1);

end
