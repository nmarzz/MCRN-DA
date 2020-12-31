function [ q ] = getDMD( Phi,p )
%Phi is the matrix of DMD modes
q = Phi(:,...
    1:min(p, size(Phi,2))...
    ); %this get us p colmuns, or the whole matrix

% DMD modes are typically complex valued. Here we replace conjugate pairs
% with real/imag parts

for k = 1:size(q,2) - 1
    
    % test if two consecutive columns are conjugates
    if norm( q(:,k) - conj(q(:,k+1)) ) < 1e-6
        
        % take cartesian components
        realpart = real(q(:,k));
        imagpart = imag(q(:,k));
        
        % replace two columns with cartesian components
        q(:,k) = realpart;
        q(:,k+1) = imagpart;
        
    end
    
end

%% modified gram schmidt to get the orthogonal basis
q = mgs(q);
