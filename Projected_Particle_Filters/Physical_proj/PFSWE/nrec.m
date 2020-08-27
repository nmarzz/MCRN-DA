function [u] = nrec(n)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    if (n == 0) | (n ==1)
        u = 1;
    else
        u = n*nrec(n-1)*nrec(n-2);
end

