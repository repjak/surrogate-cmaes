classdef GenerationsUpdaterRankDiff < GenerationsUpdater
  methods
    function err = computeErr(obj, ~, ~, ~, modelY, origY, dim, lambda, countiter, varargin)
      err = errRankMu(modelY, origY, obj.ec.cmaesState.mu);
    end

    function obj = GenerationsUpdaterRankDiff(ec, parameters)
      % constructor
      obj = obj@GenerationsUpdater(ec, parameters);
    end
  end
end
