function [q] = getpod(Ur,p)
%Ur is the reduce matrix 
[~,r]=size(Ur);
if p<=r
    q=Ur(:,1:p); %this get us p colmuns
else
    q=Ur; % this is when p>r
   
end
% size(q)