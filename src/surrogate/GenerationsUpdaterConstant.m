classdef GenerationsUpdaterConstant < GenerationsUpdater
  properties
    origGenerations
    modelGenerations
  end

  methods
    % get new values
    function [origGenerations, modelGenerations] = update(obj, ~, ~, ~, ~, ~)
      origGenerations = obj.origGenerations;
      modelGenerations = obj.modelGenerations;
    end

    function obj = GenerationsUpdaterConstant(parameters)
      obj = obj@GenerationsUpdater(parameters);

      if isstruct(parameters)
        obj.origGenerations = defopts(parameters, 'origGenerations', 1);
        obj.modelGenerations = defopts(parameters, 'modelGenerations', 1);
      elseif iscell(parameters)
        parsedParams = struct(parameters{:});
        obj.origGenerations = defopts(parsedParams, 'origGenerations', 1);
        obj.modelGenerations = defopts(parsedParams, 'modelGenerations', 1);
      else
        assert(isnumeric(parameters) && length(parameters) == 2, 'Invalid parameters for GenerationsUpdaterConstant');
        obj.origGenerations = parameters(1);
        obj.modelGenerations = parameters(2);
      end
    end
  end
end