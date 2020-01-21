function xnew = uvhtox(u,v,h,nx,ny)

u1=reshape(u,nx*ny,1);
v1=reshape(v,nx*ny,1);
h1=reshape(h,nx*ny,1);
xnew=cat(1,u1,v1,h1);
