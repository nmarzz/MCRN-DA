function explorePF(samples) % per each projection scheme (ie, data = pod, phys = dmd ...)
rng(1331);
% Test PF code, look for optimality in choosing projection dimension and other parameters
N = 40; % Dimension of Lorenz

% LHS_Call(xmin,xmean,xmax,xsd,nsample,distrib,threshold)
Ls = round(repmat(LHS_Call(500, 1000, 2000, 0 ,samples,'unif'),9,1));
alphas = repmat(LHS_Call(1e-2, 0.5, 1, 0 ,samples,'unif'),9,1);
resampThreshs = repmat(LHS_Call(1e-2, 0.5, 1, 0 ,samples,'unif'),9,1); 
tolp = round(repmat(LHS_Call(N/100, N/2, N, 0 ,samples,'unif'),9,1));
told =  round(repmat(LHS_Call(N/100, N/2, N, 0 ,samples,'unif'),9,1));
nump =  round(repmat(LHS_Call(N/100, N/2, N, 0 ,samples,'unif'),9,1));
numd =  round(repmat(LHS_Call(N/100, N/2, N, 0 ,samples,'unif'),9,1));
PhysicalProjection = [zeros(samples*3,1) ; ones(samples*3,1);2*ones(samples*3,1)];
DataProjection = repmat([zeros(samples,1) ; ones(samples,1);2*ones(samples,1)],3,1);

numruns = samples*9;
pfvals = nan(numruns,3);
for i = 1:numruns
    fprintf('\n%d/%d',i,numruns);
    
   [RMSEave_orig,RMSEave_proj,ResampPercent] = pf(N,Ls(i),alphas(i),resampThreshs(i),PhysicalProjection(i),DataProjection(i),tolp(i),told(i),nump(i),numd(i));
   pfvals(i,:,:,:) = [RMSEave_orig,RMSEave_proj,ResampPercent];
end


scatter(pfvals(:,3),pfvals(:,2)) % resampling percent vs rmse projected

save('pfexplore')
end