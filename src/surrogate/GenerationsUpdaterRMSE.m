classdef GenerationsUpdaterRMSE < GenerationsUpdater
  properties
    minModelGenerations
    maxModelGenerations
    updateRate
    rmse
    generations
    lastModelGenerations
    origGenerations
    steepness
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
    function [origGenerations, modelGenerations] = update(obj, modelY, origY, ~, ~, countiter)
      obj.generations(end+1) = countiter;

      if isempty(modelY)
        obj.rmse(end+1) = NaN;
      else      
        err = sqrt(sum((modelY - origY).^2))/length(origY);
        obj.rmse(end+1) = err;

        obj.lastModelGenerations = obj.minModelGenerations + GenerationsUpdaterRMSE.simplesig(1 - err / max(obj.rmse), obj.steepness) * (obj.maxModelGenerations - obj.minModelGenerations);

        disp(['GenerationsUpdaterRMSE: model generations: ', num2str(round(obj.lastModelGenerations)), ' [ ', repmat('+', 1, round(obj.lastModelGenerations)), repmat(' ', 1, obj.maxModelGenerations - round(obj.lastModelGenerations)), ' ]']);
      end

      origGenerations = round(obj.origGenerations);
      modelGenerations = round(obj.lastModelGenerations);
    end

    function obj = GenerationsUpdaterRMSE(parameters)
      % constructor
      obj = obj@GenerationsUpdater(parameters);
      if isstruct(parameters)
        parsedParams = parameters;
      else
        parsedParams = struct(parameters{:});
      end

      obj.minModelGenerations = defopts(parsedParams, 'minModelGenerations', 1);
      obj.maxModelGenerations = defopts(parsedParams, 'maxModelGenerations', 5);
      obj.updateRate = defopts(parsedParams, 'updateRate', 0.5);
      obj.steepness = defopts(parsedParams,'steepness', 5);
      obj.rmse = [];
      obj.origGenerations = 1;
      obj.lastModelGenerations = obj.minModelGenerations;
    end
  end
end
