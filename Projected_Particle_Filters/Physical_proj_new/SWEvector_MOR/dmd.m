function out = dmd( DataMatrix, dt, r, options )
%% DMD Dynamic Mode Decomposition (exact version)
%
% out = dmd( DataMatrix, dt, r, ... )
%
% Inputs:
%     DataMatrix (:,:) double {mustBeNumeric, mustBeReal, mustBeFinite} -- INPUT SNAPSHOT MATRIX (columns are snapshots)
%     dt (1,1) double {mustBePositive, mustBeFinite} -- INPUT SNAPSHOT TIMESTEP (columns are snapshots)
%     r (1,1) double {mustBePositive,mustBeInteger,mustBeFinite} -- REQUESTED DIMENSION OF REDUCED-ORDER MODEL
%
% Additional optional name-value pairs, specified as
% out = dmd( DataMatrix, dt, r, NAME, VALUE, NAME, VALUE,... )
%
%     out = dmd( DataMatrix, dt, r, 'rom_type', VALUE)
%         VALUE = 'lsq' - standard least-squares truncation {default}
%         VALUE = 'tlsq' - total least squares truncation (Hemati 2017)
%
%    out = dmd( DataMatrix, dt, r, 'sortbyb', VALUE)
%         VALUE = true - sort DMD modes by |b|
%         VALUE = false - sort DMD modes by average L2 norm of time
%                          coefficient {default}
%
%    out = dmd( DataMatrix, dt, r, 'removemean', VALUE)
%         VALUE = true - remove the mean before computing DMD (Hirsh 2019)
%         VALUE = false - no preprocessing {default}
%
%    out = dmd( DataMatrix, dt, r, 'step', VALUE)
%         VALUE is a row-vector of snapshot indices at which optimization is
%         performed to determine b coefficients in the expansion. 
%         If VALUE = -1, all snapshots are used.
% 
%        
%
% Outputs:
% out.meanL2norm - time-average of L2 norm of mode
% out.b - contribution of the mode to 0th timestep
% out.Phi - DMD mode vector
% out.omega - continuous time DMD eigenvalue omega = log( lambda ) / dt
% out.lambda - discrete time DMD eigenvalue lambda = exp( omega * dt )
% out.model_rank - rank of the model (= input r parameter)

arguments

    DataMatrix (:,:) double {mustBeNumeric, mustBeReal, mustBeFinite}
    dt (1,1) double {mustBePositive, mustBeFinite}
    r (1,1) double {mustBePositive,mustBeInteger,mustBeFinite}
    options.rom_type {mustBeMember(options.rom_type,{'lsq','tlsq'})} = 'lsq'
    options.removemean (1,1) logical = false
    options.sortbyb (1,1) logical = false
    options.step (1,:) double {mustBeInteger,mustBeFinite} = 1
end

if options.step == -1
    options.step = 1:size(DataMatrix,2);
end

assert( r<= min(size(DataMatrix)), ...
    "ROM order has to be smaller than rank of input")

%% Preprocessing
Nsnapshots = size(DataMatrix,2);

X1=DataMatrix(:,1:end-1); % likely a very slow part of the code
X2=DataMatrix(:,2:end);

if options.removemean % if de-meaning requested

    X1m = mean(X1, 2); % compute actual mean vectors
    X2m = mean(X2, 2);

else % just make "mean vectors" equal to zero
    X1m = zeros(size(X1, 1), 1 );
    X2m = zeros(size(X2, 1), 1 );

end

% subtract means
X1 = X1 - X1m;
X2 = X2 - X2m;

% total-DMD Hemati debias
switch(options.rom_type)
    case 'tlsq'
        disp("Total least-squares truncation to order " + r)
        Z = [X1;X2];
        [~,~,Q] = svd(Z,'econ');
        Qr = Q(:,1:r);
        X1 = X1*Qr;
        X2 = X2*Qr;
    case 'lsq'
        disp("Standard least-squares truncation to order " + r)
end


%% Core DMD algorithm

% SVD
[U, S, V] = svd(X1, 'econ');
Ur = U(:, 1:r);
Sr = S(1:r, 1:r);
Vr = V(:,1:r);

% Build Atilde
Atilde = transpose(Ur)*X2*Vr*inv(Sr);
assert( all(~isnan(Atilde),'all') && all(~isinf(Atilde),'all'), ...
    "Atilde shouldn't contain NaN or Inf - check invertibility of Sr")

assert( length(Atilde) == r, "Requested model rank not achieved")

% Compute DMD Modes
[W, D] = eig(Atilde);
Phi = X2*Vr*inv(Sr)*W;

% normalize each column by its 2-norm
Phi = Phi ./ vecnorm( Phi );

%% Compute continuous-time eigenvalues
lambda = diag(D);
omega = log(lambda)/dt;

%% compute combination coefficients
% by L2-fit to a sequence of snapshots

disp("Computing b by L2 fit at snapshots " + ...
    mat2str(options.step) + ", time = " + ...
    mat2str( (options.step-1)*dt ) );

% LHS: modes advanced to required steps
LHS = zeros([size(Phi), numel(options.step)]);

% advance mode vectors by steps and set them into layers of LHS tensor
for k = 1:numel(options.step)
    LHS(:,:,k) = Phi*D^(options.step(k)-1); 
end

% reshape LHS into a 2D matrix so that the layers are stacked on top of
% each other
LHS = permute(LHS, [1,3,2]);
LHS = reshape(LHS,[], size(Phi,2) );

% basically 
% LHS = [ Phi * D^0 ; Phi * D^1 ; ...]


% RHS: corresponding snapshots
RHS = DataMatrix(:, options.step);
RHS = reshape(RHS,[], 1);

% basically
% RHS = [Datamatrix(:,1); Datamatrix(:,2);, ...]

% not always equivalent to LHS\RHS
% see: https://www.mathworks.com/help/matlab/ref/lsqminnorm.html#mw_e9cfa6e3-ccc3-4830-802c-df496b9e452b
b = lsqminnorm(LHS, RHS);


%% Compute mean L2 contribution of each mode
T = (Nsnapshots-1)*dt;
meanL2norm = abs(b) .* sqrt( (exp(2*real(omega)*T)-1)./(2*real(omega)*T) );


%% Sorting
if options.sortbyb
    [~,idx] = sort( abs(b), 'descend' );
else
    [~,idx] = sort( meanL2norm, 'descend' );
end

out.meanL2norm = meanL2norm(idx);
out.b = b(idx);
out.Phi = Phi(:,idx);
out.omega = omega(idx);
out.lambda = lambda(idx);
out.model_rank = length(Atilde);
end
