function [A,R] = mgs(A)
%Performs modified Gram-Schmidt to find Q,R
%such that A=QR.
%A is over written with orthgonal matrix Q

[m, n] = size(A);
R = zeros(n,n);
for j= 1:n
    R(j,j)=norm(A(:,j));
    A(:,j)=A(:,j)/R(j,j);
    R(j,j+1:n)=A(:,j)'*A(:,j+1:n);
    A(:,j+1:n)=A(:,j+1:n)-A(:,j)*R(j,j+1:n);
end
