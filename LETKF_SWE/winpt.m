function [window,point]=winpt(k,irad,wlen,n)
%This works with non-periodic BCs
      if k-irad >= 1 & k+irad <=n
         window = circshift([1:n]', irad-k+1); %will use center point
         point = irad+1;
      elseif k+irad <=n
         window = [1:n]'; %will use "kth" point
         point = k;
      else
         window = [n-2*irad:n]'; %will use kth point
         point=wlen+k-n;
      end
      %point = irad+1;

