function x=Hty(y,inth,N,L)
% Transpose of H
x=zeros(N,L);
x(1:inth:end,:) = y;
end