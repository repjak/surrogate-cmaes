classdef GenerationsUpdaterKL < GenerationsUpdater
  properties
    cmConstants
    lastMu
    lastLambda
    historyKL
    % the historical maximum can be overestimated by a discount factor at normalization to avoid too pessimistic behaviour
    normDiscountFactor
    % a number of iterations in which the discount factor linearly diminishes to 1
    normDiscountIterations
  end

  methods
    function err = computeErr(obj, arx, arxvalid, arz, modelY, origY, ~, ~, countiter, varargin)
      [xmean1, C1, sigma1] = cmaesUpdate(arx, arxvalid, arz, modelY, obj.ec.cmaesState, obj.cmConstants);
      [xmean2, C2, sigma2] = cmaesUpdate(arx, arxvalid, arz, origY, obj.ec.cmaesState, obj.cmConstants);
      newKL = mvnKL(xmean1, sigma1*C1, xmean2, sigma2*C2);

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

    function [origGenerations, modelGenerations] = update(obj, arx, arxvalid, arz, modelY, origY, dim, mu, lambda, countiter, varargin)
      if mu ~= obj.lastMu || lambda ~= obj.lastLambda
        % cmaes internal parameters (from optimalized cmaes code)
        % TODO: perhaps initialize these as part of cmOptions
        chiN = dim^0.5*(1-1/(4*dim)+1/(21*dim^2));
        weights = log(max(mu, lambda/2) + 1/2)-log(1:mu)';
        weights = weights/sum(weights); % normalize recombination weights
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

      obj.lastMu = mu;
      obj.lastLambda = lambda;

      [origGenerations, modelGenerations] = update@GenerationsUpdater(obj, arx, arxvalid, arz,  modelY, origY, dim, mu, lambda, countiter, varargin{:});
    end

    function obj = GenerationsUpdaterKL(ec, parameters)
      % constructor
      obj = obj@GenerationsUpdater(ec, parameters);
	    obj.normDiscountFactor = defopts(obj.parsedParams, 'geneECAdaptive_normDiscountFactor', 2);
	    obj.normDiscountIterations = defopts(obj.parsedParams, 'geneECAdaptive_normDiscountIterations', 5);
      obj.cmConstants = struct();
      obj.lastMu = NaN;
      obj.lastLambda = NaN;
    end
  end
end
