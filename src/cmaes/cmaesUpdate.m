function [xmean, C, sigma] = cmaesUpdate(arx, arfitness, cmaesState, cmaesConstants)
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
  diagD = cmaesState.diagD;
  ps = cmaesState.ps;
  pc = cmaesState.pc;
  sigma = cmaesState.sigma;
  counteval = cmaesState.countiter;
  N = cmaesState.dim;

  invsqrtC = B * diag(diagD.^(-1)) * B';

  % Sort by fitness and compute weighted mean into xmean
  [~, arindex] = sort(arfitness);  % minimization
  xold = xmean;
  xmean = arx(:,arindex(1:mu)) * weights;  % recombination, new mean value

  % Cumulation: Update evolution paths
  ps = (1-cs) * ps ...
        + sqrt(cs*(2-cs)*mueff) * invsqrtC * (xmean-xold) / sigma;
  hsig = sum(ps.^2)/(1-(1-cs)^(2*counteval/lambda))/N < 2 + 4/(N+1);
  pc = (1-cc) * pc ...
        + hsig * sqrt(cc*(2-cc)*mueff) * (xmean-xold) / sigma;

  % Adapt covariance matrix C
  artmp = (1/sigma) * (arx(:,arindex(1:mu)) - repmat(xold,1,mu));  % mu difference vectors
  C = (1-c1-cmu) * C ...                   % regard old matrix
       + c1 * (pc * pc' ...                % plus rank one update
               + (1-hsig) * cc*(2-cc) * C) ... % minor correction if hsig==0
       + cmu * artmp * diag(weights) * artmp'; % plus rank mu update

  % Adapt step size sigma
  sigma = sigma * exp((cs/damps)*(norm(ps)/chiN - 1));
end
