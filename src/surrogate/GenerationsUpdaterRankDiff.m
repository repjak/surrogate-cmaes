classdef GenerationsUpdaterRankDiff < GenerationsUpdater
  methods
    function err = computeError(obj, modelY, origY, dim, lambda, countiter, varargin)
      [~, sort1] = sort(modelY);
      ranking2   = ranking(origY);
      err = errRankMuOnly(ranking2(sort1), obj.ec.cmaesState.mu);
    end

    function obj = GenerationsUpdaterRankDiff(ec, parameters)
      % constructor
      obj = obj@GenerationsUpdater(parameters);
    end
  end
end
