function [j,k] = mattovec(i,nx,ny)

%i = j + (k-1)*nx;
k = ceil(i/nx);
j = i - (k-1)*nx;
