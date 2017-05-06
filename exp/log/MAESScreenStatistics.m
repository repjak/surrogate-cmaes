classdef MAESScreenStatistics < Observer
%SCREENSTATISTICS -- print statistics from ModelAssistedEC on screen
  properties
    verbosity
  end

  methods
    function obj = MAESScreenStatistics(params)
      obj@Observer();
      verbosity = defopts(params, 'verbose', 5);
    end

    function notify(obj, ec, varargin)
      % get the interesting data and process them
      if (mod(ec.cmaesState.countiter, 10) == 1)
        fprintf('######## iter /evals | D_fopt. | rmse | Kendall        | M  nData | sigma^2. |\n');
      end
      model = '.';
      nTrainData = 0;
      if (~isempty(ec.model) && ec.model.isTrained() ...
          && ec.model.trainGeneration == ec.cmaesState.countiter)
        model = '+'; nTrainData = ec.model.getTrainsetSize(); end
      outputValues1 = [...
          ec.cmaesState.countiter, ec.counteval, ...
          ec.stats.fmin - ec.surrogateOpts.fopt, ...
          ec.stats.rmse, ...
          ec.stats.kendall ];
      outputValues2 = [...
          nTrainData, ...
          ec.cmaesState.sigma^2 ];
      outputValues1(isnan(outputValues1)) = 0.0;
      outputValues2(isnan(outputValues2)) = 0.0;
      %         #####  iter /evals(or,p) | Dopt |rmseR | rnkR | rnk2 |rnkVal * | Mo nD nDiR |sigm^2| aErr |smooEr| orRat| aGain|
      fprintf('=[MAES]= %4d /%5d | %.1e | %.2f | %.2f %s | %s %2d | %.2e |\n', ... %  %.2f | %.2f | %.2f | %.2f |\n', ...
          outputValues1(:), decorateKendall(ec.stats.kendall), model, outputValues2(:) ...
          );
    end
  end
end
