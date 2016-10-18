classdef GenerationsUpdaterRMSEKendall < GenerationsUpdater
  properties
    rmse
    kendall
    alpha
  end

  methods
    function err = computeErr(obj, ~, ~, ~, modelY, origY, ~, ~, countiter, varargin)
      newRMSE = sqrt(sum((modelY - origY).^2))/length(origY);
      newKendall = corr(modelY, origY, 'type', 'Kendall');

      obj.rmse(countiter) = newRMSE;
      obj.kendall(countiter) = newKendall;

      maxRMSE = max(obj.rmse);

      % combine the RMSE and Kendall ranking coefficient
      err = (1 - obj.alpha) * (newRMSE / maxRMSE) + obj.alpha * (1 - newKendall) / 2;
    end

    function obj = GenerationsUpdaterRMSEKendall(ec, parameters)
      % constructor
      obj = obj@GenerationsUpdater(ec, parameters);
      % weight of Kendall coeficient over RMSE, from (0, 1)
      obj.alpha = defopts(obj.parsedParams, 'geneECAdaptive_weightKendallRMSE', 1);
    end
  end
end
