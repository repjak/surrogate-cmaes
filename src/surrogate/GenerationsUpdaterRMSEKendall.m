classdef GenerationsUpdaterRMSEKendall < GenerationsUpdater
  properties
    parsedParams
    minModelGenerations
    maxModelGenerations
    errThreshold
    updateRate
    rmse
    kendall
    err
    generations
    lastModelGenerations
    origGenerations
    steepness
    alpha
    transferFun
  end

  methods (Static)
    function y = simplesig(x, k)
      % a simple "normalized" sigmoid function obtaining values 0 and 1 in
      % points 0 and 1, respectively and having the middle point in 0.5
      % @k -- a steepness factor, k > 0
      r = 1 + 1/k;
      y = r * (x - 0.5) ./ (r - 1 + abs(2 * (x - 0.5))) + 0.5;
    end
  end

  methods
    function [newOrigGenerations, newModelGenerations] = update(obj, modelY, origY, ~, ~, countiter)
      obj.generations(end+1) = countiter;

      if isempty(modelY)
        obj.rmse(end+1) = NaN;
        obj.kendall(end+1) = NaN;
        obj.err(end+1) = NaN;
      else
        newRMSE = sqrt(sum((modelY - origY).^2))/length(origY);
        newKendall = corr(modelY, origY, 'type', 'Kendall');

        if ~isempty(obj.rmse) && ~all(isnan(obj.rmse))
          maxRMSE = max(obj.rmse);
        else
          maxRMSE = newRMSE;
        end

        % combine the RMSE and Kendall ranking coefficient
        newErr = (1 - obj.alpha) * (newRMSE / maxRMSE) + obj.alpha * (1 - newKendall) / 2;

        % exponential smoothing
        if ~isempty(obj.err) && ~all(isnan(obj.err))
          errValid = obj.err(~isnan(obj.err));
          lastErr = errValid(end);
        else
          lastErr = 0.5;
        end

        newErr = (1 - obj.updateRate) * lastErr + obj.updateRate * newErr;

        obj.rmse(end+1) = newRMSE;
        obj.kendall(end+1) = newKendall;
        obj.err(end+1) = newErr;

        % normalize into an interval below an error threshold and apply a
        % transfer function
        errNormalized = obj.transferFun(max(1 - newErr / obj.errThreshold, 0));

        % min-max scaling and rounding
        obj.lastModelGenerations = round(obj.minModelGenerations + errNormalized * (obj.maxModelGenerations - obj.minModelGenerations));

        disp(['GenerationsUpdaterRMSE: model generations: ', num2str(round(obj.lastModelGenerations)), ' [ ', repmat('+', 1, round(obj.lastModelGenerations)), repmat(' ', 1, obj.maxModelGenerations - round(obj.lastModelGenerations)), ' ]']);
      end

      newOrigGenerations = round(obj.origGenerations);
      newModelGenerations = round(obj.lastModelGenerations);
    end

    function obj = GenerationsUpdaterRMSEKendall(parameters)
      % constructor
      obj = obj@GenerationsUpdater(parameters);
      if ~isstruct(parameters)
        obj.parsedParams = struct(parameters{:});
      else
        obj.parsedParams = parameters;
      end
      obj.minModelGenerations = defopts(obj.parsedParams, 'minModelGenerations', 1);
      obj.maxModelGenerations = defopts(obj.parsedParams, 'maxModelGenerations', 5);
      obj.updateRate = defopts(obj.parsedParams, 'updateRate', 0.5);
      obj.errThreshold = defopts(obj.parsedParams, 'errThreshold', 0.45);
      obj.steepness = defopts(obj.parsedParams, 'steepness', 5);
      obj.alpha = defopts(obj.parsedParams, 'alpha', 0.2);
      obj.transferFun = defopts(obj.parsedParams, 'transferFun', '@(x) x');
      obj.rmse = [];
      obj.err = [];
      obj.origGenerations = 1;
      obj.lastModelGenerations = obj.minModelGenerations;
    end
  end
end