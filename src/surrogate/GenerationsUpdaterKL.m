classdef GenerationsUpdaterKL < GenerationsUpdater
  properties
    cmConstants
    historyKL
    % the historical maximum can be overestimated by a discount factor at normalization to avoid too pessimistic behaviour
    normDiscountFactor
    % a number of iterations in which the discount factor linearly diminishes to 1
    normDiscountIterations
  end

  methods (Static)
    function divKL = mvnKL(m1, C1, m2, C2)
      % compute Kullback-Leibler divergence of two multivariate
      % distributions given by means 'm1' and 'm2', respectively
      % and positive semidefinite covariance matrices 'C1' and 'C2',
      % respectively

      assert(length(m1) == length(m2) && all(size(C1) == size(C2)), 'Dimensions don''t match');

      L1 = chol(C1, 'lower');
      logdetC1 = GenerationsUpdaterKL.cov_logdet(L1);

      L2 = chol(C2, 'lower');
      logdetC2 = GenerationsUpdaterKL.cov_logdet(L2);
      invC2 = GenerationsUpdaterKL.cov_inv(L2);

      dim = length(m1);
      trprod = invC2(:)'*C1(:);
      mdiff = m2 - m1;

      divKL = 0.5 * (trprod + mdiff' * invC2 * mdiff - dim + (logdetC2 - logdetC1));
      divKL = max(0, divKL);
    end

    function cov_logdet = cov_logdet(L)
      % compute log(det(C)) for a positive semidefinite matrix given by
      % a lowertriangular matrix from a LU factorization of C
      cov_logdet = 2 * sum(log(diag(L)));
    end

    function cov_inv = cov_inv(L)
      % compute covariance matrix inverse from its lower triangular LU-factor
      invL = L\eye(size(L, 1));
      cov_inv = invL' * invL;
    end
  end

  methods
    function err = computeErr(obj, arxvalid, modelY, origY, ~, ~, countiter, varargin)
      [xmean1, C1, sigma1] = cmaesUpdate(obj, arxvalid, modelY, obj.ec.cmaesState);
      [xmean2, C2, sigma2] = cmaesUpdate(obj, arxvalid, origY, obj.ec.cmaesState);
      newKL = obj.mvnKL(xmean1, sigma1*C1, xmean2, sigma2*C2);

      obj.historyKL(end+1:(countiter-1)) = NaN;
      obj.historyKL(countiter) = newKL;
      maxKL = max(obj.historyKL);

      if maxKL == 0
        err = 0;
        return;
      end

      if (obj.normDiscountIterations == 0 || obj.normDiscountFactor < 1)
        discount = 1;
      else
        it = length(obj.historyKL(~(obj.historyKL == 0) & ~isnan(obj.historyKL))) - 1;
        discount = max(obj.normDiscountFactor - it * (obj.normDiscountFactor - 1) / obj.normDiscountIterations, 1);
      end
      err = newKL / (discount *  maxKL);
    end

    function [xmean, C, sigma] = cmaesUpdate(obj, arx, arfitness, cmaesState)
      % one-step update of CMA-ES distribution
      % based on http://www.lri.fr/~hansen/purecmaes.m

      weights = obj.cmConstants.weights;
      mueff = obj.cmConstants.mueff;
      cc = obj.cmConstants.cc;
      cs = obj.cmConstants.cs;
      damps = obj.cmConstants.damps;
      chiN = obj.cmConstants.chiN;
      c1 = obj.cmConstants.c1;
      cmu = obj.cmConstants.cmu;

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
      dim = cmaesState.dim;
      N = dim;

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

    function [origGenerations, modelGenerations] = update(obj, arxvalid, modelY, origY, dim, mu, lambda, countiter, varargin)
      if isempty(fields(obj.cmConstants))
        % initialization (from optimalized cmaes)
        % TODO: perhaps initialize these as part of cmOptions
        chiN = dim^0.5*(1-1/(4*dim)+1/(21*dim^2));
        weights = log(max(mu, lambda/2) + 1/2)-log(1:mu)';
        mueff = sum(weights)^2/sum(weights.^2); % variance-effective size of mu
        cc = (4 + mueff/dim) / (dim+4 + 2*mueff/dim); % cumulation constant for pc'
        cs = (mueff+2)/(dim+mueff+3);
        damps = 1 + 2*max(0,sqrt((mueff-1)/(dim+1))-1) + cs;
        c1 = 2 / ((dim+1.3)^2+mueff);  % learning rate for rank-one update'
        cmu = 2 * (mueff-2+1/mueff) / ((dim+2)^2+mueff); % learning rate for rank-mu update'

        obj.cmConstants = struct( ...
          'chiN', chiN, ...
          'weights', weights, ...
          'mueff', mueff, ...
          'cc', cc, ...
          'cs', cs, ...
          'damps', damps, ...
          'c1', c1, ...
          'cmu', cmu ...
        );
      end

      [origGenerations, modelGenerations] = update@GenerationsUpdater(obj, arxvalid, modelY, origY, dim, mu, lambda, countiter, varargin{:});
    end

    function obj = GenerationsUpdaterKL(ec, parameters)
      % constructor
      obj = obj@GenerationsUpdater(ec, parameters);
	    obj.normDiscountFactor = defopts(obj.parsedParams, 'geneECAdaptive_normDiscountFactor', 2);
	    obj.normDiscountIterations = defopts(obj.parsedParams, 'geneECAdaptive_normDiscountIterations', 5);
      obj.cmConstants = struct();
    end
  end
end

