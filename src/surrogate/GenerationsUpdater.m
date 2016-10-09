classdef (Abstract) GenerationsUpdater < handle
  % In a generation-based evolution control, a regression model is used for
  % certain number of generations, which we can call the model's lifelength.
  %
  % In generations when new data are available and a new model is trained,
  % method 'update' estimates (adapts) the new model's lifelength based on
  % the previous model's error.
  %
  % The output is given as a pair [newOrigGenerations, newModelGenerations],
  % where newModelGenerations indicates for how long the new model should
  % be used, while newOrigGenerations is recommended number of
  % fitness-evaluated generations that follow (usually constant).
  properties
    parsedParams
    ec
    maxModelGenerations
    minModelGenerations
    startModelGenerations
    lastModelGenerations
    origGenerations
    updateRate
    lowErrThreshold
    highErrThreshold
    aggregateType
    weights
    errSmoothingFactor
    errSmoothed
    transferFun
    gain
    plotDebug
    historyErr
    historyAggErr
    historyModelGenerations
    lastUpdateIteration
    fh
  end

  methods (Abstract)
    err = computeErr(obj, arxvalid, modelY, origY, dim, lambda, countiter, varargin);
  end

  methods (Static)
    function y = simplesig(x, k)
      % a simple "normalized" sigmoid function obtaining values 0 and 1 in
      % points 0 and 1, respectively and having the middle point in 0.5
      % @k -- a steepness factor, k > 0
      r = 1 + 1/k;
      y = r * (x - 0.5) ./ (r - 1 + abs(2 * (x - 0.5))) + 0.5;
    end

    function Kmn = covFcn(Xm, Xn, theta)
      % A covariance function defined as sum of SEard, constant
      % shift and a white noise:
      %
      % Kmn(x_i, x_j) = sf2 * exp(-(x_i - x_q)'*diag(1./ell)*(x_i - x_j)/2)
      %                 + sc2
      %                 + sn2*delta(x_i, x_j)
      %
      % where ell = (ell_1^2,...ell_D^2) is vector of ARD parameters, D is
      % dimension of input space, sf2 is signal variance, sc2 is constant shift
      % and sn2 is added white noise.
      %
      % Theta is a vector of positive numbers:
      % theta = [ log(ell_1)
      %           log(ell_2)
      %            .
      %           log(ell_D)
      %           log(sqrt(sf2))
      %           log(sqrt(sc2))
      %           log(sqrt(sn2)) ]
      %
      assert(size(Xm, 2) == size(Xn, 2), 'Dimensions must agree');
      dim = size(Xm, 2);

      % tolerance on squared distance of vectors when they're considered equal
      tol = 1e-9;

      assert(all(size(theta) == [dim+3 1]), ['theta must be a row vector of ''' dim + 3 ''' hyperparameters']);

      ell = exp(theta(1:dim));
      sf2 = exp(2*theta(dim+1));
      sc2 = exp(2*theta(dim+2));
      sn2 = exp(2*theta(dim+3));

      Kmn = sq_dist(diag(1./ell) * Xm', diag(1./ell) * Xn');
      Kmn = sf2 * exp(-Kmn/2);
      Kmn = bsxfun(@plus, Kmn, sc2);

      delta = bsxfun(@lt, sq_dist(Xm', Xn'), tol^2);
      Kmn(delta) = Kmn(delta) + sn2;
    end
  end

  methods
    function [origGenerations, modelGenerations] = update(obj, arxvalid, modelY, origY, dim, mu, lambda, countiter, varargin)
      % a template method to estimate new model lifelength
      if (nargin >= 9), obj.ec = varargin{1}{:}; end
      obj.historyErr((obj.lastUpdateIteration+1):(countiter-1)) = NaN;
      obj.historyModelGenerations((obj.lastUpdateIteration+1):(countiter-1)) = obj.lastModelGenerations;

      % compute current error and add it to history
      if (isempty(modelY) || (max(modelY) - min(modelY)) == 0 ...
          || (max(origY) - min(origY)) == 0)
        err = NaN;
      else
        err = computeErr(obj, arxvalid, modelY, origY, dim, lambda, countiter);
      end
      obj.historyErr(countiter) = err;

      % aggregate historical errors
      aggErr = obj.aggregateWithHistory();

      % the error is truncated by a threshold interval
      % a transfer function is applied
      % the gain is scaled into an admissible interval
      if (~isnan(aggErr))
        obj.gain = min(max(0, aggErr - obj.lowErrThreshold), (obj.highErrThreshold - obj.lowErrThreshold)) / (obj.highErrThreshold - obj.lowErrThreshold);
        obj.gain = obj.transferFun(1-obj.gain);
        newModelGenerations = obj.minModelGenerations + obj.gain * (obj.maxModelGenerations - obj.minModelGenerations);
      else
        obj.gain = NaN;
        newModelGenerations = obj.startModelGenerations;
      end
      % Debug:
      fprintf('err = %.2f ;  gain = %.2f ; output = %.2f\n', err, obj.gain, newModelGenerations);

      % obj.lastModelGenerations -- the last used number of model generations (from the last updated orig generation)
      if (countiter > 1)
        obj.lastModelGenerations = obj.historyModelGenerations(countiter-1);
      else
        obj.lastModelGenerations = obj.startModelGenerations;
      end

      % the final model's lifelength is exponentially updated:   l = (1-a) * l_old  +  a * l_new
      modelGenerations = round((1-obj.updateRate) * obj.lastModelGenerations + obj.updateRate * newModelGenerations);
      modelGenerations = min(max(modelGenerations, obj.minModelGenerations), obj.maxModelGenerations);

      obj.historyModelGenerations(countiter) = modelGenerations;
      obj.lastModelGenerations = modelGenerations;
      obj.lastUpdateIteration = countiter;

      if obj.plotDebug
        fprintf('New model generations=%0.2f based on err trend=%0.2f\n', modelGenerations, newModelGenerations);
        obj.historyAggErr(countiter) = aggErr;
      end

      origGenerations = obj.origGenerations;
    end

    function value = aggregateWithHistory(obj)
      % aggregate last criterion values into one value
      %
      % This implementation:
      % - ignores NaN values in history (less values are used then)
      %
      % (a) takes median of the last length(weights) values
      %     or
      % (b) takes weighted sum of the last length(weights) values

      % take at most nHistory last values
      nHistory = min(length(obj.historyErr), length(obj.weights));
      values = obj.historyErr((end-nHistory+1):end);
      % identify NaN's
      bValues  = ~isnan(values);
      % return with NaN if no valid values in history
      if (~any(bValues))
        value = NaN;
        return;
      end

      switch lower(obj.aggregateType)
      case 'median'
        value = median(values(bValues));
      case 'weightedsum'
        % take adequate weights and re-norm them to sum to 1
        localWeights = obj.weights((end-nHistory+1):end);
        localWeights = localWeights(bValues) ./ sum(localWeights(bValues));
        % return the weighted sum
        value = sum(localWeights .* values(bValues));
      case 'lastvalid'
        % take the last valid (non-NaN) value
        nonNanValues = values(bValues);
        value = nonNanValues(end);
      case 'last'
        % take the last value (even if NaN)
        value = obj.historyErr(end);
      case 'expsmoothing'
        if length(bValues) < 5
          obj.errSmoothed = median(values(bValues));
          value = obj.errSmoothed;
        else
          lastbValue = bValues(end);
          obj.errSmoothed = (1-obj.errSmoothingFactor)*obj.errSmoothed + obj.errSmoothingFactor*values(lastbValue);
          value = obj.errSmoothed;
        end
      otherwise
        error('GenerationsUpdater: aggregateType ''%s'' not implemented.', obj.aggregateType);
      end
    end

    function obj = GenerationsUpdater(ec, parameters)
      % constructor
      if ~isstruct(parameters)
        obj.parsedParams = struct(parameters{:});
      else
        obj.parsedParams = parameters;
      end
      % parameter 'ec' is a reference to the EvolutionControl
      obj.ec = ec;
      % lower and upper bounds of model's lifelength, maxModelGenerations > minModelGenerations >= 0
      obj.minModelGenerations = defopts(obj.parsedParams, 'geneECAdaptive_minModelGenerations', 1);
      obj.maxModelGenerations = defopts(obj.parsedParams, 'geneECAdaptive_maxModelGenerations', 5);
      obj.startModelGenerations = defopts(obj.parsedParams, 'geneECAdaptive_startModelGenerations', round(obj.maxModelGenerations - obj.minModelGenerations)/2);
      obj.lastModelGenerations = obj.startModelGenerations;
      % weights for weighted sum of historical error values
      obj.weights = defopts(obj.parsedParams, 'geneECAdaptive_weights', exp((1:4)/2) / sum(exp((1:4)/2)));
      % smoothing factor for historical error values
      obj.errSmoothingFactor = defopts(obj.parsedParams, 'geneECAdaptive_errSmoothingFactor', 0.5);
      % exponential smoothing of output values
      obj.updateRate = defopts(obj.parsedParams, 'geneECAdaptive_updateRate', 0.9);
      % type of aggregation of historical values of RankDiff errors
      obj.aggregateType = defopts(obj.parsedParams, 'geneECAdaptive_aggregateType', 'weightedSum');
      % a transfer function applied to error, such as a sigmoid function
      obj.transferFun = defopts(obj.parsedParams, 'geneECAdaptive_transferFun', @(x) x);

      % lowest and highest error which affect gain util it saturates to 0 or 1
      obj.lowErrThreshold  = defopts(obj.parsedParams, 'geneECAdaptive_lowErrThreshold', 0.1);
      obj.highErrThreshold = defopts(obj.parsedParams, 'geneECAdaptive_highErrThreshold', 0.5);

      obj.historyErr = [];
      obj.historyAggErr = [];
      obj.origGenerations = 1;
      obj.lastUpdateIteration = 0;
      obj.plotDebug = 0;
    end
  end
end
