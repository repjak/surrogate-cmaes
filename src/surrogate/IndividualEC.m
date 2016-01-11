classdef IndividualEC
  properties
    model
  end
  
  methods
    function obj = IndividualEC()
    % constructor
      obj.model = [];
    end
    
    function [fitness_raw, arx, arxvalid, arz, counteval, lambda, archive, surrogateStats] = runGeneration(obj, cmaesState, surrogateOpts, archive, varargin)
      % Run one generation of individual evolution control
      
      fitness_raw = [];
      arx = [];
      arxvalid = [];
      arz = [];
      
      xmean = cmaesState.xmean;
      sigma = cmaesState.sigma;
      lambda = cmaesState.lambda;
      BD = cmaesState.BD;
      diagD = cmaesState.diagD;
      fitfun_handle = cmaesState.fitfun_handle;
      countiter = cmaesState.countiter;
      counteval = surrogateOpts.sampleOpts.counteval;
      
      obj.model = ModelFactory.createModel(surrogateOpts.modelType, surrogateOpts.modelOpts, xmean');

      if (isempty(obj.model))
        % model could not be created :( use the standard CMA-ES
        return;
      end
      
      minTrainSize = obj.model.getNTrainData();

      nArchivePoints = myeval(surrogateOpts.evoControlTrainNArchivePoints);
      [xTrain, yTrain] = archive.getDataNearPoint(nArchivePoints, ...
          xmean', surrogateOpts.evoControlTrainRange, sigma, BD);
      
      [fitness_raw, arx, arxvalid, arz, archive, counteval, xTrain, yTrain] = ...
        presample(minTrainSize, surrogateOpts, cmaesState, archive, xTrain, yTrain, varargin{:});

      % train the model
      % TODO: omit the unnecessary variables xmean, sigma and BD
      % as they are already in cmaesState  
      obj.model = obj.model.train(xTrain, yTrain, xmean', countiter, sigma, BD, cmaesState);

      missingTrainSize = max(minTrainSize - size(xTrain, 1), 0);
      nLambdaRest = lambda - missingTrainSize;
      
      if (~obj.model.isTrained())
        return
      end
      
      if any(strcmpi(obj.model.predictionType, {'poi', 'ei'}))
        bestImprovement = 0;
        % sample 'gamma' populations of size 'nLambdaRest'
        for sampleNumber = 1:surrogateOpts.evoControlIndividualExtension
          [xExtend, xExtendValid, zExtend] = ...
              sampleCmaesNoFitness(xmean, sigma, nLambdaRest, BD, diagD, surrogateOpts.sampleOpts);
          % TODO: criterion for choosing the best sample
          actualImprovement = mean(obj.model.getModelOutput(xExtend'));
          % choose sample with higher improvement factor (PoI, EI)
          if actualImprovement > bestImprovement || sampleNumber == 1
            xToReeval = xExtend;
            xToReevalValid = xExtendValid;
            zToReeval = zExtend;
            bestImprovement = actualImprovement;
          end
        end

      else
        % sample the enlarged population of size 'gamma * nLambdaRest'
        extendSize = ceil(surrogateOpts.evoControlIndividualExtension ...
            * nLambdaRest);
        [xExtend, xExtendValid, zExtend] = ...
            sampleCmaesNoFitness(xmean, sigma, extendSize, BD, diagD, surrogateOpts.sampleOpts);
        % calculate the model prediction for the extended population
        yExtend = obj.model.getModelOutput(xExtend');

        nBest = min(ceil(lambda*surrogateOpts.evoControlBestFromExtension), nLambdaRest - 1);
        nCluster = nLambdaRest - nBest;
        [xToReeval, xToReevalValid, zToReeval] = ...
            SurrogateSelector.choosePointsToReevaluate(...
            xExtend, xExtendValid, zExtend, yExtend, nBest, nCluster);
      end
      
      % original-evaluate the chosen points
      [yNew, xNew, xNewValid, zNew, counteval] = ...
          sampleCmaesOnlyFitness(xToReeval, xToReevalValid, zToReeval, xmean, sigma, nLambdaRest, BD, diagD, fitfun_handle, surrogateOpts.sampleOpts, varargin{:});
      surrogateOpts.sampleOpts.counteval = counteval;
      fprintf('counteval: %d\n', counteval)
      % update the Archive
      archive = archive.save(xNewValid', yNew', countiter);
      % the obj.models' dataset will be supplemented with this
      % new points during the next training using all the xTrain
      % calculate the models' precision
      yPredict = obj.model.predict(xNewValid');
      kendall = corr(yPredict, yNew', 'type', 'Kendall');
      rmse = sqrt(sum((yPredict' - yNew).^2))/length(yNew);
      fprintf('  model-gener.: %d preSamples, reevaluated %d pts, test RMSE = %f, Kendl. corr = %f.\n', missingTrainSize, nLambdaRest, rmse, kendall);
      surrogateStats = [rmse kendall];
      
      % save the resulting re-evaluated population as the returning parameters
      fitness_raw = [fitness_raw yNew];
      arx = [arx xNew];
      arxvalid = [arxvalid xNewValid];
      arz = [arz zNew];

    end
    
  end
  
end