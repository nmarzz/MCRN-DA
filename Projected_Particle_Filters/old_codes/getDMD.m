function [ q ] = getDMD( Phi,p )
%Phi is the matrix of DMD modes
[~,r] = size(Phi);
if p<=r
    q = Phi(:,1:p); %this get us p colmuns
    q=mgs(q);
else
    q = Phi; % this is when p>r
    q=mgs(q);
end

