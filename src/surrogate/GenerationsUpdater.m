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
    updateRateDown
    defaultErr
    lowErrThreshold
    highErrThreshold
    transferFun
    gain
    historyErr
    historySmoothedErr
    historyModelGenerations
    lastUpdateIteration
  end

  methods (Abstract)
    err = computeErr(obj, arx, arxvalid, arz, modelY, origY, dim, lambda, countiter, varargin);
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
    function [origGenerations, modelGenerations] = update(obj, arx, arxvalid, arz, modelY, origY, dim, mu, lambda, countiter, varargin)
      % a template method to estimate new model lifelength
      if (nargin >= 9), obj.ec = varargin{:}; end
      obj.historyErr((obj.lastUpdateIteration+1):(countiter-1)) = NaN;
      obj.historySmoothedErr((obj.lastUpdateIteration+1):(countiter-1)) = NaN;
      obj.historyModelGenerations((obj.lastUpdateIteration+1):(countiter-1)) = obj.lastModelGenerations;

      % compute current error and add it to history
      if (isempty(modelY) || (max(modelY) - min(modelY)) == 0 ...
          || (max(origY) - min(origY)) == 0)
        err = NaN;
      else
        err = computeErr(obj, arx, arxvalid, arz, modelY, origY, dim, lambda, countiter);
        if isnan(err)
          fprintf('GenerationsUpdater.update(): current error is NaN');
        end
      end
      obj.historyErr(countiter) = err;

      % calculate exponentialy smoothed value of the error
      %      e_{t} = (1-a) * e_{t-1}  +  a * err
      % and add it to the history
      lastIdx = find(~isnan(obj.historySmoothedErr(1:max(countiter-1,0))), 1, 'last');
      if (isempty(lastIdx))
        % there's no valid smoothed error value in history, use 'defaultErr'
        lastSmoothedErr = obj.defaultErr;
      else
        % take the last non-NaN smoothed error value
        lastSmoothedErr = obj.historySmoothedErr(lastIdx);
      end
      if (~isnan(err))
        % we have got a new current error ==> use EWA smooth
        if (err > lastSmoothedErr)
          % err is higher than last smoothed
          smoothedErr = (1-obj.updateRate) * lastSmoothedErr + obj.updateRate * err;
        else
          % err is lower than last smoothed
          smoothedErr = (1-obj.updateRateDown) * lastSmoothedErr + obj.updateRateDown * err;
        end

        obj.historySmoothedErr(countiter) = smoothedErr;
      else
        % we do not have a new error value, just use the last one
        smoothedErr = lastSmoothedErr;
        obj.historySmoothedErr(countiter) = NaN;
      end

      % the error is truncated by a threshold interval
      % a transfer function is applied
      % the gain is scaled into an admissible interval
      obj.gain = min(max(0, smoothedErr - obj.lowErrThreshold), (obj.highErrThreshold - obj.lowErrThreshold)) / (obj.highErrThreshold - obj.lowErrThreshold);
      obj.gain = obj.transferFun(1-obj.gain);
      modelGenerations = round(obj.minModelGenerations + obj.gain * (obj.maxModelGenerations - obj.minModelGenerations));

      % Debug:
      fprintf('err = %.2f ;  gain = %.2f ; output = %.2f\n', err, obj.gain, modelGenerations);

      obj.historyModelGenerations(countiter) = modelGenerations;
      obj.lastModelGenerations = modelGenerations;
      obj.lastUpdateIteration = countiter;

      origGenerations = obj.origGenerations;
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
      % exponential smoothing of error values
      obj.updateRate = defopts(obj.parsedParams, 'geneECAdaptive_updateRate', 0.40);
      % for negative updates, use the same rate as positive rate as default;
      % this default is used also if [] is supplied in experiment definition
      % for 'updateRateDown'
      obj.updateRateDown = defopts(obj.parsedParams, 'geneECAdaptive_updateRateDown', obj.updateRate);
      % a transfer function applied to error, such as a sigmoid function
      obj.transferFun = defopts(obj.parsedParams, 'geneECAdaptive_transferFun', '@(x) x');
      if ischar(obj.transferFun)
        obj.transferFun = evalin('caller', obj.transferFun);
      end

      % lowest and highest error which affect gain util it saturates to 0 or 1
      obj.lowErrThreshold  = defopts(obj.parsedParams, 'geneECAdaptive_lowErrThreshold', 0.10);
      obj.highErrThreshold = defopts(obj.parsedParams, 'geneECAdaptive_highErrThreshold', 0.50);

      obj.defaultErr = defopts(obj.parsedParams, 'geneECAdaptive_defaultErr', 0.20);

      obj.historyErr = [];
      obj.historySmoothedErr = [];
      obj.origGenerations = 1;
      obj.lastUpdateIteration = 0;
    end
  end
end
