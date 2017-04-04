classdef DoubleTrainedEC < EvolutionControl & Observable
%
% TODO:
% [ ] remove updaterParams and use DTAdaptive_* parameters instead
% [ ] rename 'restrictedParam' to 'origRatio'
% [ ] in choosePointsForReevaluation() consider lowering 'mu' according to the proportion size(xExtend,2) / obj.cmaesState.lambda
%
  properties 
    model
    pop
    cmaesState
    counteval

    origRatioUpdater
    restrictedParam
    useDoubleTraining
    maxDoubleTrainIterations
    retrainedModel
    stats
    usedUpdaterState            % Updater's state variables updated in the last generation
    archive
    nPresampledPoints
    surrogateOpts
    newModel
    modelArchive
    modelArchiveGenerations
    modelArchiveLength
    acceptedModelAge            % how many generations old model is still OK (0 == only current)
    modelAge                    % age of model in the number of generations (0 == current model)
    oldModelAgeForStatistics    % age of model for gathering statistics of old models
    isTrainSuccess
    origPointsRoundFcn % function computing number of original-evaluated points from origRatio
    nBestPoints                 % the number of points with the best predicted f-value to take every generation
    usedBestPoints              % how many best-predicted points was really orig-evaluated
    validationGenerationPeriod  % the number of generations between "validation generations" + 1, validation generation is used when 'mod(g, validationGenerationPeriod) == 0'; see validationPopSize
    validationPopSize           % the minimal number of points to be orig-evaluated in validation generation
  end

  methods 
    function obj = DoubleTrainedEC(surrogateOpts, varargin)
    % constructor
      obj@Observable();
      obj.model = [];
      obj.pop = [];
      obj.surrogateOpts = surrogateOpts;

      % DTS parameters
      obj.restrictedParam = defopts(surrogateOpts, 'evoControlRestrictedParam', 0.1);
      obj.useDoubleTraining = defopts(surrogateOpts, 'evoControlUseDoubleTraining', true);
      obj.maxDoubleTrainIterations = defopts(surrogateOpts, 'evoControlMaxDoubleTrainIterations', Inf);

      % other initializations:
      obj.acceptedModelAge = defopts(surrogateOpts, 'evoControlAcceptedModelAge', 2);
      obj.origPointsRoundFcn = str2func(defopts(surrogateOpts, 'evoControlOrigPointsRoundFcn', 'ceil'));

      % Adaptive DTS parameters
      surrogateOpts.updaterType = defopts(surrogateOpts, 'updaterType', 'none');
      surrogateOpts.updaterParams = defopts(surrogateOpts, 'updaterParams', {});
      obj.origRatioUpdater = OrigRatioUpdaterFactory.createUpdater(obj, surrogateOpts);

      % Preselection parameters
      obj.nBestPoints = defopts(surrogateOpts, 'evoControlNBestPoints', 0);
      obj.validationGenerationPeriod = defopts(surrogateOpts, 'evoControlValidationGenerationPeriod', 1);
      obj.validationPopSize = defopts(surrogateOpts, 'evoControlValidationPopSize', 0);

      % Model Archive fixed settings and initialization
      obj.oldModelAgeForStatistics = [3:5];
      obj.modelArchiveLength = 5;
      obj.modelArchive = cell(1, obj.modelArchiveLength);
      obj.modelArchiveGenerations = nan(1, obj.modelArchiveLength);
      obj.modelAge = 0;
      obj.isTrainSuccess = false;

      % statistics
      obj.stats = struct( ...
          'fmin', NaN, ...              % minimal original fitness in population
          'rmseReeval', NaN, ...        % RMSE of the re-evaluated point(s)
          'kendallReeval', NaN, ...     % Kendall's corr. of the re-evaluated point(s)
          'rankErrReeval', NaN, ...     % rank error of popul. with re-evaluated point(s)
          'rankErr2Models', NaN, ...    % rank error between prediction of two models
          'rmseValid', NaN, ...         % RMSE of the (2nd) model on the validation set
          'kendallValid', NaN, ...      % Kendall of the (2nd) model on the validation set
          'rankErrValid', NaN, ...      % rank error between true fitn. and model pred. on the validation set
          'rmseOldModel', NaN, ...      % RMSE of old model on future orig points
          'kendallOldModel', NaN, ...   % Kendall of old model on future orig points
          'normKendallOldModel', NaN, ... % Kendall of old model normed to [0,1]
          'ageOldModel', NaN, ...       % age of the old model which used for statistics
          'nDataOldModel', 0, ...       % the number of data points from archive for old model statistics
          'lastUsedOrigRatio', NaN, ... % restricted param which was used (last) in the last generation
          'adaptErr', NaN, ...     % last measured rankDiff during update()
          'adaptGain', NaN, ...         % gain of original ratio (to be converted via min/max)
          'adaptSmoothedErr', NaN ...   % smoothed error value used before fed into transfer function
          );
      obj.usedUpdaterState = struct( ...
          'gain', NaN, ...
          'err', NaN, ...
          'smoothedErr', NaN ...
          );
    end

    function [obj, fitness_raw, arx, arxvalid, arz, counteval, lambda, archive, surrogateStats, origEvaled] = runGeneration(obj, cmaesState, surrogateOpts, sampleOpts, archive, counteval, varargin)
    % Run one generation of double trained evolution control

      % initialization
      lambda = cmaesState.lambda;       % this is needed due to myeval
      dim    = cmaesState.dim;          % this is needed due to myeval
      obj.cmaesState = cmaesState;
      obj.archive    = archive;
      obj.counteval  = counteval;
      obj.retrainedModel = [];
      obj.stats.nDataInRange = NaN;
      obj.modelAge = 0;
      obj.isTrainSuccess = false;
      obj.usedBestPoints = 0;

      % prepare the final population to be returned to CMA-ES
      obj.pop = Population(lambda, dim);
      obj.restrictedParam = obj.origRatioUpdater.update([], [], dim, lambda, obj.cmaesState.countiter, obj);

      obj.newModel = ModelFactory.createModel(obj.surrogateOpts.modelType, obj.surrogateOpts.modelOpts, obj.cmaesState.xmean');

      if (isempty(obj.newModel))
        [obj, ok] = obj.tryOldModel();
        if (~ok)
          % model could not be created nor older is usable :(. Use the standard CMA-ES.
          [obj, fitness_raw, arx, arxvalid, arz, counteval, surrogateStats, origEvaled] ...
              = obj.finalizeGeneration(sampleOpts, varargin);
          return;
        end
      end

      % if a new model is used, find appropriate training set and train it:
      if (obj.modelAge == 0)

        minTrainSize = obj.newModel.getNTrainData();

        nArchivePoints = myeval(obj.surrogateOpts.evoControlTrainNArchivePoints);
        [xTrain, yTrain, nData] = obj.archive.getDataNearPoint(nArchivePoints, ...
            obj.cmaesState.xmean', obj.surrogateOpts.evoControlTrainRange, ...
            obj.cmaesState.sigma, obj.cmaesState.BD);
        obj.stats.nDataInRange = nData;

        % Do pre-sample
        [ok, y, arx, x, arz, ~, obj.counteval, xTrain, yTrain] = ...
            presample(minTrainSize, obj.cmaesState, obj.surrogateOpts, sampleOpts, ...
            obj.archive, obj.counteval, xTrain, yTrain, varargin{:});
        obj.nPresampledPoints = size(x, 2);
        phase = 0;        % pre-sampled points
        obj.pop = obj.pop.addPoints(x, y, arx, arz, obj.nPresampledPoints, phase);

        if (~ok)
          % not enough data for training model
          [obj, ok] = obj.tryOldModel();
          if (~ok)
            [obj, fitness_raw, arx, arxvalid, arz, counteval, surrogateStats, origEvaled] ...
                = obj.finalizeGeneration(sampleOpts, varargin);
            return;
          end
        end

        % (first) model training
        obj.newModel = obj.newModel.train(xTrain, yTrain, obj.cmaesState, sampleOpts);
        if (obj.newModel.isTrained())
          obj = obj.updateModelArchive(obj.newModel, obj.modelAge);
        else
          [obj, ok] = obj.tryOldModel();
          if (~ok)
            % model cannot be trained :( -- return with orig-evaluated population
            [obj, fitness_raw, arx, arxvalid, arz, counteval, surrogateStats, origEvaled] ...
                = obj.finalizeGeneration(sampleOpts, varargin);
            return;
          end
        end

      end  % if (obj.modelAge == 0)

      obj.model = obj.newModel;
      nLambdaRest = lambda - obj.nPresampledPoints;

      % Validation Generation -- raise the number of orig. evaluated points
      % once in several (opts.validationGenerationPeriod) generations
      if (obj.validationPopSize > 0 && mod(obj.cmaesState.countiter, obj.validationGenerationPeriod) == 0)
        obj.restrictedParam = max(min(1, obj.validationPopSize/nLambdaRest), obj.restrictedParam);
      end

      % the number of points to orig-evaluate
      nPoints = obj.origPointsRoundFcn(nLambdaRest * obj.restrictedParam);

      % Preselection: orig-evaluate the best predicted point(s)
      % out of 50*lambda sampled points (if obj.nBestPoints > 0).
      % Saves the really used number of preselected points into obj.usedBestPoints.
      [obj, yBestOrig, xBest, xBestValid, zBest] = obj.preselection( ...
          obj.nBestPoints, nPoints, sampleOpts, varargin{:});
      nPoints = nPoints - obj.usedBestPoints;

      % sample new points -- the rest of population
      [xExtend, xExtendValid, zExtend] = ...
          sampleCmaesNoFitness(obj.cmaesState.sigma, nLambdaRest - obj.usedBestPoints, obj.cmaesState, sampleOpts);

      % merge the preselected best point(s) with the sampled rest of population
      xExtend      = [xBest,      xExtend];
      xExtendValid = [xBestValid, xExtendValid];
      zExtend      = [zBest,      zExtend];
      isEvaled = [true(1,obj.usedBestPoints), false(1, nLambdaRest-obj.usedBestPoints)];
      yOrig    = [yBestOrig, NaN(1, nLambdaRest-obj.usedBestPoints)];

      % get the model's prediction
      [modelOutput, yExtendModel] = obj.model.getModelOutput(xExtendValid');

      doubleTrainIteration = 0;
      notEverythingEvaluated = true;

      % main DTS re-evaluating cycle
      while (notEverythingEvaluated)

        doubleTrainIteration = doubleTrainIteration + 1;

        % save statistics about these just used values
        obj.stats.lastUsedOrigRatio = obj.restrictedParam;
        obj.stats.adaptGain = obj.usedUpdaterState.gain;
        obj.stats.adaptErr = obj.usedUpdaterState.err;
        obj.stats.adaptSmoothedErr = obj.usedUpdaterState.smoothedErr;

        if (nPoints > 0)

          % choose point(s) for re-evaluation
          reevalID = false(1, nLambdaRest);
          reevalID(~isEvaled) = obj.choosePointsForReevaluation(nPoints, ...
              xExtend(:, ~isEvaled), modelOutput(~isEvaled), yExtendModel(~isEvaled));
          xToReeval = xExtendValid(:, reevalID);
          nToReeval = sum(reevalID);

          % original-evaluate the chosen points
          [yNew, xNew, xNewValid, zNew, obj.counteval] = ...
              sampleCmaesOnlyFitness(xExtend(:, reevalID), xToReeval, zExtend(:, reevalID), ...
              obj.cmaesState.sigma, nToReeval, obj.counteval, obj.cmaesState, sampleOpts, ...
              varargin{:});
          xExtendValid(:, reevalID) = xNewValid;
          xExtend(:, reevalID) = xNew;
          zExtend(:, reevalID) = zNew;
          yOrig(reevalID) = yNew;
          isEvaled = isEvaled | reevalID;
          % Debug:
          % fprintf('counteval: %d\n', obj.counteval)

          phase = 1;        % re-evaluated points
          obj.pop = obj.pop.addPoints(xNewValid, yNew, xNew, zNew, nToReeval, phase);

          % update the Archive
          obj.archive.save(xNewValid', yNew', obj.cmaesState.countiter);

        else
          yNew = []; xNew = []; xNewValid = []; zNew = [];
        end % if (nPoints > 0)

        if (sum(isEvaled) > 0)

          % re-train the model again with the new original-evaluated points
          xTrain = [xTrain; xNewValid'];
          yTrain = [yTrain; yNew'];
          obj.retrainedModel = obj.model.train(xTrain, yTrain, obj.cmaesState, sampleOpts);

          if (obj.useDoubleTraining && obj.retrainedModel.isTrained())
            % if internal CMA-ES restart just happend, create a new OrigRatioUpdater
            if (~isempty(obj.origRatioUpdater.lastUpdateGeneration) ...
                && obj.origRatioUpdater.lastUpdateGeneration > obj.cmaesState.countiter)
              obj.origRatioUpdater = OrigRatioUpdaterFactory.createUpdater(obj, obj.surrogateOpts);
              obj = obj.updateModelArchive(obj.retrainedModel, obj.modelAge);
            end

            % origRatio adaptivity (ratio will be used for the next iteration
            % of the while cycle or for the next generation of CMA-ES)
            yFirstModel  = obj.model.predict(xExtendValid');
            yExtendModel = obj.retrainedModel.predict(xExtendValid');
            obj.restrictedParam = obj.origRatioUpdater.update(...
                yFirstModel', yExtendModel', dim, lambda, obj.cmaesState.countiter, obj);

            % save the statistics about the new Updater's state
            obj.usedUpdaterState.err = obj.origRatioUpdater.historyErr(obj.cmaesState.countiter);
            obj.usedUpdaterState.gain = obj.origRatioUpdater.gain;
            obj.usedUpdaterState.smoothedErr = obj.origRatioUpdater.historySmoothedErr(obj.cmaesState.countiter);
          end

        end % if (sum(isEvaled) > 0)

        % determine the number of points to re-evalute for the next
        % iteration
        nPoints = obj.origPointsRoundFcn(nLambdaRest * obj.restrictedParam) - sum(isEvaled);

        notEverythingEvaluated = (doubleTrainIteration < obj.maxDoubleTrainIterations) ...
            && (floor(lambda * obj.restrictedParam) > sum(isEvaled));
      end % while (notEverythingEvaluated)

      if (~all(isEvaled))
        phase = 2;      % model-evaluated rest of the population
        obj.pop = obj.pop.addPoints(xExtendValid(:, ~isEvaled), yExtendModel(~isEvaled), ...
            xExtend(:, ~isEvaled), zExtend(:, ~isEvaled), 0, phase);
      end

      assert(obj.pop.nPoints == lambda, 'There are not yet all lambda points prepared, but they should be!');

      % sort the returned solutions (the best to be first)
      [obj.pop, ~] = obj.pop.sort;

      % save the resulting re-evaluated population as the returning parameters
      [obj, fitness_raw, arx, arxvalid, arz, counteval, surrogateStats, origEvaled] ...
          = obj.finalizeGeneration(sampleOpts, varargin);
    end


    function obj = updateModelArchive(obj, newModel, modelAge)
      % update the modelArchive with the current new model
      countiter = obj.cmaesState.countiter;

      if (obj.modelArchiveGenerations(1) < (countiter - modelAge))
        % there's an old model in the first position ==> shift old
        % models to the history
        obj.modelArchive(2:end) = obj.modelArchive(1:(end-1));
        obj.modelArchiveGenerations(2:end) = obj.modelArchiveGenerations(1:(end-1));
        % clear the first position
        obj.modelArchive{1} = [];
        obj.modelArchiveGenerations(1) = NaN;
        obj.isTrainSuccess = true;
      end

      if (newModel.isTrained())
        % the obj.newModel should be usable and newer than what we have
        % in the first position, so save it there
        obj.modelArchive{1} = newModel;
        obj.modelArchiveGenerations(1) = (countiter - modelAge);
      end
    end


    function [obj, ok] = tryOldModel(obj)
      % try to load an appropriate old model
      % if successful, load it as the newModel
      [oldModel, age] = obj.getOldModel(0:obj.acceptedModelAge);
      if (~isempty(oldModel))
        obj.newModel = oldModel;
        obj.modelAge = age;
        ok = true;
      else
        ok = false;
      end
    end


    function [oldModel, modelAge] = getOldModel(obj, generationDiffs)
      % return an old model from modelArchive from the generations
      %   countiter - [generationDiffs]
      % Returns the first such model found, or [] if none is found
      %
      % Examples:
      % m = dec.getOldModel(3:5)        % returns the youngest model 3--5 generations old
      % m = dec.getOldModel(0)          % returns the model from current generation

      countiter = obj.cmaesState.countiter;
      modelGenerations = countiter - generationDiffs;
      modelGenerations = modelGenerations(modelGenerations > 0);

      oldModel = []; modelAge = NaN;
      if (isempty(modelGenerations))
        % no sensible generations given, return []
        return;
      end

      % go through modelArchive and look whether the models (in
      % historical order) are not from the specified generations
      for i = 1:obj.modelArchiveLength
        ind = find(modelGenerations == obj.modelArchiveGenerations(i), 1, 'first');
        if (ind)
          oldModel = obj.modelArchive{i};
          modelAge = (countiter - obj.modelArchiveGenerations(i));
          return;
        end
      end
      % no model from specified generations found, return []
    end


    function [obj, fitness_raw, arx, arxvalid, arz, counteval, surrogateStats, origEvaled] = finalizeGeneration(obj, sampleOpts, varargin)
      % fill the rest of the population with original evaluations
      [obj, fitness_raw, arx, arxvalid, arz, counteval] = ...
          obj.fillPopWithOrigFitness(sampleOpts, varargin);
      origEvaled = obj.pop.origEvaled;

      % calculate statistics
      obj.stats.rmseReeval     = NaN; % RMSE of the re-evaluated point(s)
      obj.stats.kendallReeval  = NaN; % Kendall's corr. of the re-evaluated point(s)
      obj.stats.rankErrReeval  = NaN; % rank error of popul. with re-evaluated point(s)
      obj.stats.rankErr2Models = NaN; % rank error between prediction of two models
      obj.stats.rmseValid      = NaN; % RMSE of the (2nd) model on the validation set
      obj.stats.kendallValid   = NaN; % Kendall of the (2nd) model on the validation set
      obj.stats.rankErrValid   = NaN; % rank error between true fitness and model pred.
      obj.stats.rmseOldModel   = NaN; % RMSE of old model on future orig points
      obj.stats.kendallOldModel = NaN; % Kendall of old model on future orig points
      obj.stats.normKendallOldModel = NaN; % Kendall transformed to [0,1] error-range
      obj.stats.ageOldModel    = NaN;
      obj.stats.nDataOldModel  = 0;
      obj.stats.fmin = min(obj.pop.y(1,(obj.pop.origEvaled == true)));

      [obj.stats.rmseOldModel, obj.stats.kendallOldModel, ...
          obj.stats.normKendallOldModel, obj.stats.ageOldModel, ...
          obj.stats.nDataOldModel] = obj.oldModelStatistics();

      % model-related statistics
      if (~isempty(obj.model) && obj.model.isTrained() ....
          && ~all(obj.pop.origEvaled))
        % predict the population by the first model
        yModel1 = obj.model.predict(obj.pop.x');

        if (any(obj.pop.origEvaled))
          % calculate RMSE, Kendall's coeff. and ranking error
          % between the original fitness and the first model's values
          % of the re-evaluated point(s), i.e. (phase == 1)
          [obj.stats.rmseReeval, obj.stats.kendallReeval, obj.stats.rankErrReeval] ...
              = obj.reevalStatistics(yModel1);

          % get ranking error between the first and the second model
          % (if the second model is trained)
          obj.stats.rankErr2Models = obj.retrainStatistics(yModel1);
        end

        % independent validation set statistics
        [obj.stats.rmseValid, obj.stats.kendallValid, obj.stats.rankErrValid, lastModel] ...
            = obj.validationStatistics(sampleOpts);

        % shift the f-values:
        %   if the model predictions are better than the best original value
        %   in the model's dataset, shift ALL (!) function values
        %   Note: - all values have to be shifted in order to preserve predicted
        %           ordering of values
        %         - small constant is added because of the rounding errors
        %           when numbers of different orders of magnitude are summed
        fminDataset = min(lastModel.dataset.y);
        fminModel = obj.pop.getMinModeled;
        diff = max(fminDataset - fminModel, 0);
        obj.pop = obj.pop.shiftY(1.000001*diff);
        fitness_raw = obj.pop.y;
      end

      obj.notify_observers();

      % for backwards compatibility
      [minstd minstdidx] = min(obj.cmaesState.sigma*sqrt(obj.cmaesState.diagC));
      [maxstd maxstdidx] = max(obj.cmaesState.sigma*sqrt(obj.cmaesState.diagC));
      surrogateStats = [obj.stats.rmseValid, obj.stats.kendallValid, obj.cmaesState.sigma, ...
        max(obj.cmaesState.diagD)/min(obj.cmaesState.diagD), minstd, maxstd, ...
        obj.stats.rankErr2Models];
    end


    function [obj, yNew, xNew, xNewValid, zNew, counteval] = fillPopWithOrigFitness(obj, sampleOpts, varargin)
      %Fill the rest of the current population 'pop' (of class Population) with 
      % the original-evaluated individuals
      nToEval = obj.cmaesState.lambda - sum(obj.pop.nPoints);
      if (nToEval > 0)
        [y, arx, x, arz, obj.counteval] = sampleCmaes(obj.cmaesState, sampleOpts, ...
            nToEval, obj.counteval, varargin{:});
        obj.archive.save(x', y', obj.cmaesState.countiter);

        phase = 3;      % original-evaluated rest of the population
        obj.pop = obj.pop.addPoints(x, y, arx, arz, nToEval, phase);
      end
      obj.pop.sort();
      yNew = obj.pop.y;
      xNewValid = obj.pop.x;
      xNew = obj.pop.arx;
      zNew = obj.pop.arz;
      counteval = obj.counteval;
    end

    function reevalID = choosePointsForReevaluation(obj, nPoints, xExtend, modelOutput, yExtendModel)
    % choose points with low confidence to reevaluate
    %
    % returns:
    %   reevalID  --  bool vector which points from xExtend to reevaluate
      if any(strcmpi(obj.model.predictionType, {'sd2', 'poi', 'ei'}))
        % higher criterion is better (sd2, poi, ei)
        [~, pointID] = sort(modelOutput, 'descend');

      elseif (strcmpi(obj.model.predictionType, 'expectedrank'))
        MIN_POINTS_FOR_EXPECTED_RANK = 3;

        ok = true;
        if (~isempty(obj.retrainedModel) && obj.retrainedModel.isTrained())
          lastModel = obj.retrainedModel;
        elseif (~isempty(obj.model) && obj.model.isTrained())
          lastModel = obj.model;
        else
          warning('No valid model for calculating expectedRankDiff(). Using "sd2" criterion.');
          ok = false;
        end
        if (size(xExtend, 2) < MIN_POINTS_FOR_EXPECTED_RANK)
          fprintf(2, 'expectedRankDiff(): #pop=%d < %d: using sd2 criterion\n', size(xExtend, 2), MIN_POINTS_FOR_EXPECTED_RANK);
          ok = false;
        end
        % TODO: consider lowering 'mu' according to the proportion
        %       size(xExtend,2) / obj.cmaesState.lambda
        if (ok)
          [pointID, errs] = expectedRankDiff(lastModel, xExtend, obj.cmaesState.mu);
          if (~ sum(errs >= eps) > (size(xExtend,2)/2))
            warning('exptectedRankDiff() returned more than lambda/2 points with zero expected rankDiff error. Using "sd2" criterion.');
            ok = false;
          end
        end
        if (~ok)
          [~, sd2] = lastModel.predict(xExtend');
          [~, pointID] = sort(sd2, 'descend');
        end
        % Debug:
        % y_r = ranking(y_m);
        % fprintf('  Expected permutation of sorted f-values: %s\n', num2str(y_r(pointID)'));

      else
        % lower criterion is better (fvalues, lcb, fpoi, fei)
        [~, pointID] = sort(yExtendModel, 'ascend');
      end

      nLambdaRest = size(xExtend, 2);
      reevalID = false(1, nLambdaRest);
      assert(obj.origRatioUpdater.getLastRatio() >= 0 && obj.origRatioUpdater.getLastRatio() <= 1, 'origRatio out of bounds [0,1]');
      reevalID(pointID(1:nPoints)) = true;
    end


    function [rmse, kendall, rankErr] = reevalStatistics(obj, yModel1)
      % calculate RMSE and possibly Kendall's coeff. of the re-evaluated point(s)
      % (phase == 1)
      phase1 = (obj.pop.phase == 1);
      if (~any(phase1))
        rmse = NaN; kendall = NaN; rankErr = NaN;
        return;
      end
      yReeval = obj.pop.y(1,phase1);
      yReevalModel = yModel1(phase1);
      rmse = sqrt(sum((yReevalModel' - yReeval).^2))/length(yReeval);
      kendall = corr(yReevalModel, yReeval', 'type', 'Kendall');

      yModel2 = yModel1;
      yModel2(phase1) = yReeval;
      rankErr = errRankMu(yModel1, yModel2, obj.cmaesState.mu);

      % Debug:
      % fprintf('  model: %d pSmpls, reeval %d pts, RMSE= %.2e, Kendl= %.2f, rankErr= %.3f\n', ...
      %     obj.nPresampledPoints, length(yReeval), rmse, kendall, rankErr);
    end


    function [rmse, kCorr, normKendall, age, nData] = oldModelStatistics(obj)
      % return statistics of a several-generations-old model measured
      % on new data from following generations (generations after the
      % model was trained)
      rmse = NaN; kCorr = NaN; normKendall = NaN; age = NaN; nData = 0;

      [oldModel, age] = obj.getOldModel(obj.oldModelAgeForStatistics);
      if (~isempty(oldModel))
        countiter = obj.cmaesState.countiter;
        oldModelGeneration = oldModel.trainGeneration;
        testGenerations = [(oldModel.trainGeneration+1):countiter];
        [xOrig, yOrig] = obj.archive.getDataFromGenerations(testGenerations);
        nData = length(yOrig);
        if (nData > 0)
          yPredict = oldModel.predict(xOrig);
          rmse = sqrt(sum((yPredict - yOrig).^2))/length(yPredict);
          if (nData >= 2)
            kCorr = corr(yPredict, yOrig, 'type', 'Kendall');
            normKendall = (-kCorr + 1) / 2;
          end
        end
      end
    end

    function rankErr = retrainStatistics(obj, yModel1)
    % get ranking error between the first and the second model
    % (if the second model is trained)
      if (~isempty(obj.retrainedModel) && obj.retrainedModel.isTrained())
        yModel1AfterPresample = yModel1(obj.pop.phase ~= 0);
        xAfterPresample = obj.pop.x(:,obj.pop.phase ~= 0);
        yModel2AfterPresample = obj.retrainedModel.predict(xAfterPresample');
        rankErr = errRankMu(yModel1AfterPresample, yModel2AfterPresample, obj.cmaesState.mu);
        % Debug:
        % fprintf('  2 models rank error: %.3f                  %s\n', rankErr, decorateKendall(1-rankErr*2));
      else
        % ranking error between two models cannot be calculated :(
        rankErr = NaN;
      end
    end


    function [rmse, kendall, errRank, lastModel] = validationStatistics(obj, sampleOpts)
    % generate independent validation population of points,
    % evaluate them with the original BBOB function (if we know it)
    % and calculate statistics on this independent set
    %
    % lastModel -- return the last valid model (either first or second)

      rmse = NaN; kendall = NaN; errRank = NaN;
      % do we have access to the original BBOB fitness?
      if (~ isfield(obj.surrogateOpts.modelOpts, 'bbob_func'))
        return;
      end
      if (~isempty(obj.retrainedModel) && obj.retrainedModel.isTrained())
        % statistics according to the retrained model
        lastModel = obj.retrainedModel;
      elseif (~isempty(obj.model) && obj.model.isTrained())
        % statistics according to the first model
        lastModel = obj.model;
      else
        % we do not have any valid model
        return;
      end

      [~, xValidTest, ~] = ...
          sampleCmaesNoFitness(obj.cmaesState.sigma, obj.cmaesState.lambda, obj.cmaesState, sampleOpts);
      preciseModel = ModelFactory.createModel('bbob', obj.surrogateOpts.modelOpts, obj.cmaesState.xmean');
      yTest = preciseModel.predict(xValidTest');
      yPredict = lastModel.predict(xValidTest');

      kendall = corr(yPredict, yTest, 'type', 'Kendall');
      rmse = sqrt(sum((yPredict - yTest).^2))/length(yPredict);

      errRank = errRankMu(yTest, yPredict, obj.cmaesState.mu);

      % Debug:
      % fprintf('  test RMSE= %.2e, Kendall= %.3f, rankErr= %.3f %s\n', ...
      %     rmse, kendall, errRank, decorateKendall(kendall));
    end


    function [obj, yBestOrig, xBest, xBestValid, zBest] = preselection(obj, nBestPoints, nPoints, sampleOpts, varargin)
    % Preselection: orig-evaluate the best predicted point(s)
    % out of 50*lambda sampled points

      % determine the number of presampled 'best' points
      if (nPoints <= 2 || length(nBestPoints) == 1)
        % if there _is not_ enough points to orig-evalute, use the preselection with
        % different probability...
        obj.usedBestPoints = getProbNumber(nBestPoints(1));
      else
        % ... that if there _is_ enough points
        obj.usedBestPoints = getProbNumber(nBestPoints(2));
      end

      if (nPoints >= 1 && obj.usedBestPoints >= 1)
        % The preselection should really happen.
        % select the points:
        obj.usedBestPoints = min(obj.usedBestPoints, nPoints);
        [~, xBest, xBestValid, zBest] = preselect(obj.usedBestPoints, obj.cmaesState, obj.model, sampleOpts, 50);
        % use the original function for them
        [yBestOrig,  xBest, xBestValid, zBest, obj.counteval] = ...
            sampleCmaesOnlyFitness(xBest, xBestValid, zBest, ...
            obj.cmaesState.sigma, obj.usedBestPoints, obj.counteval, obj.cmaesState, sampleOpts, ...
            varargin{:});
        % and save them
        phase = 1;        % this is the first model-evaluated point
        obj.pop = obj.pop.addPoints(xBestValid, yBestOrig, xBest, zBest, obj.usedBestPoints, phase);
        obj.archive.save(xBestValid', yBestOrig', obj.cmaesState.countiter);
      else
        % no preselection
        yBestOrig = []; xBest = []; xBestValid = []; zBest = [];
        obj.usedBestPoints = 0;
      end
    end
  end

end

function res=myeval(s)
  if ischar(s)
    res = evalin('caller', s);
  else
    res = s;
  end
end

function probNum = getProbNumber(exactNumber)
% Calculates randomized non-negative integer as follows:
%   probNum = floor(exactNumber) + eps, 
% where eps is 0 or 1. Probability that eps is 1 is equal to the remainder: 
%   P[eps = 1] = exactNumber - floor(exactNumber).
  fracN = exactNumber - floor(exactNumber);
  plus = 0;
  if (fracN > 0)
    plus = (rand() < fracN);
  end
  probNum = floor(exactNumber) + plus;
end
