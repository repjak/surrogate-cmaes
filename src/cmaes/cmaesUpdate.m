function [xmean, C, sigma] = cmaesUpdate(arx, ~, arz, arfitness, cmaesState, cmaesConstants)
  % One-step update of CMA-ES distribution and the search step.
  % Based on http://www.lri.fr/~hansen/purecmaes.m

  weights = cmaesConstants.weights;
  mueff = cmaesConstants.mueff;
  cc = cmaesConstants.cc;
  cs = cmaesConstants.cs;
  damps = cmaesConstants.damps;
  chiN = cmaesConstants.chiN;
  c1 = cmaesConstants.c1;
  cmu = cmaesConstants.cmu;

  xmean = cmaesState.xmean;
  mu = cmaesState.mu;
  lambda = cmaesState.lambda;
  C = cmaesState.C;
  B = cmaesState.B;
  ps = cmaesState.ps;
  pc = cmaesState.pc;
  sigma = cmaesState.sigma;
  counteval = cmaesState.countiter;
  N = cmaesState.dim;

  % Sort by fitness and compute weighted mean into xmean
  [~, arindex] = sort(arfitness);  % minimization
  xold = xmean;
  xmean = arx(:,arindex(1:mu)) * weights;  % recombination, new mean value
  zmean = arz(:,arindex(1:mu)) * weights;  % == D^-1*B'*(xmean-xold)/sigma

  % Cumulation: Update evolution paths
  ps = (1-cs) * ps ...
        + sqrt(cs*(2-cs)*mueff) * (B*zmean);
  hsig = sum(ps.^2)/(1-(1-cs)^(2*counteval/lambda))/N < 2 + 4/(N+1);
  pc = (1-cc) * pc ...
        + hsig * (sqrt(cc*(2-cc)*mueff)/sigma) * (xmean-exxold);

  % Adapt covariance matrix C
  artmp = (arx(:,arindex(1:mu)) - repmat(xold,1,mu)) / sigma;  % mu difference vectors

  C = (1-c1-cmu+(1-hsig)*c1*cc*(2-cc)) * C ... % regard old matrix
      + c1 * pc*pc' ...     % plus rank one update
      + cmu ...             % plus rank mu update
      * artmp * (repmat(weights,1,N) .* artmp');

  % Adapt step size sigma
  sigma = sigma * exp(min(1, (cs/damps)*(norm(ps)/chiN - 1)));
end
