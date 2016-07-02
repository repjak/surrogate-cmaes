classdef (Abstract) GenerationsUpdater < handle
  properties (Abstract)
  end
  
  methods (Abstract)
    % get new values
    [origGenerations, modelGenerations] = update(obj, modelY, origY, dim, lambda, countiter);
  end
  
  methods
    function obj = GenerationsUpdater(parameters)
      % constructor
    end
  end
end