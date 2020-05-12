function [u,v,h] = swe_vectortodata(xval, nx, ny)

x = xval;
u1=x(1:nx*ny);
v1=x(nx*ny+1:2*nx*ny);
h1=x(2*nx*ny+1:3*nx*ny);
u=reshape(u1,nx,ny);
v=reshape(v1,nx,ny);
h=reshape(h1,nx,ny);
