classdef GenerationsUpdaterRMSE < GenerationsUpdater
  properties
    parsedParams
    minModelGenerations
    maxModelGenerations
    updateRate
    rmse
    kendall
    generations
    lastModelGenerations
    origGenerations
    steepness
    alpha
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
      else      
        lastRMSE = sqrt(sum((modelY - origY).^2))/length(origY);
        lastKendall = corr(modelY, origY, 'type', 'Kendall');

        obj.rmse(end+1) = lastRMSE;
        obj.kendall(end+1) = lastKendall;

        err = (1 - obj.alpha) * (lastRMSE / max(obj.rmse)) + obj.alpha * (1 - lastKendall) / 2;

        newGenerations = obj.minModelGenerations + GenerationsUpdaterRMSE.simplesig(1 - log10(1+9*err), obj.steepness) * (obj.maxModelGenerations - obj.minModelGenerations);
        newGenerations = (1 - obj.updateRate) * obj.lastModelGenerations + obj.updateRate * newGenerations;

        obj.lastModelGenerations = min(obj.maxModelGenerations, max(obj.minModelGenerations, newGenerations))

        disp(['GenerationsUpdaterRMSE: model generations: ', num2str(round(obj.lastModelGenerations)), ' [ ', repmat('+', 1, round(obj.lastModelGenerations)), repmat(' ', 1, obj.maxModelGenerations - round(obj.lastModelGenerations)), ' ]']);
      end

      newOrigGenerations = round(obj.origGenerations);
      newModelGenerations = round(obj.lastModelGenerations);
    end

    function obj = GenerationsUpdaterRMSE(parameters)
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
      obj.steepness = defopts(obj.parsedParams, 'steepness', 5);
      obj.alpha = defopts(obj.parsedParams, 'alpha', 0.2);
      obj.rmse = [];
      obj.origGenerations = 1;
      obj.lastModelGenerations = obj.minModelGenerations;
    end
  end
end
