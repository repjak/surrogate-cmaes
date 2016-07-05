classdef (Abstract) GenerationsUpdater < handle
  properties (Abstract)
    parsedParams
  end
  
  methods (Abstract)
    % get new values
    [newOrigGenerations, newModelGenerations] = update(obj, modelY, origY, dim, lambda, countiter);
  end
  
  methods
    function obj = GenerationsUpdater(parameters)
      % constructor
    end
  end
end