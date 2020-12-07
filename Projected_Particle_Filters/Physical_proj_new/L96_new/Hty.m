function x=Hty(y,inth,N,L,minidx,maxidx)
% Transpose of H
x=zeros(N,L);
x(minidx:inth:maxidx,:) = y;
end