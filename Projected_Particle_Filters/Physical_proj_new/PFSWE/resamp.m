function [weights,X_update,numRS,Neff] = resamp(P,X,fracval)

    Neff = 1/(sum(P.^2)); % Test effective particle size
    P = P./sum(P); % Ensure particle weights are normalised

    [N,L]=size(X);

    if Neff < fracval*L
        X_update = zeros(N,L);
        CSP = cumsum(P);
        for i = 1:L
            X_update(:,i) = X(:,find(rand <= CSP,1));
        end
        weights = ones(L,1)*1./L;
        numRS=1;
    else
        weights = P;
        X_update = X;
        numRS=0;
    end
end
