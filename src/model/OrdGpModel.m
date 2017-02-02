classdef OrdGpModel < Model
  properties    % derived from abstract class "Model"
    dim                  % dimension of the input space X (determined from x_mean)
    trainGeneration = -1 % # of the generation when the model was built
    trainMean            % mean of the generation when the model was trained
    trainSigma           % sigma of the generation when the model was trained
    trainBD              % BD of the generation when the model was trained
    dataset              % .X and .y
    useShift             % whether use shift during generationUpdate()
    shiftMean            % vector of the shift in the X-space
    shiftY               % shift in the f-space
    predictionType       % type of prediction (f-values, PoI, EI)
    transformCoordinates % whether use transformation in the X-space
    stateVariables       % variables needed for sampling new points as CMA-ES do
    
    % OrdGpModel specific fields
    trainLikelihood       % negative logarithm of likelihood reached by training process
  end

  properties (Access = protected)
    stdY                  % standard deviation of Y in training set, for output normalization
    covFcn                % an identifier of covariance function
    hyp                   % a struct of hyperparameters with fields 'cov', 'ordreg', 'lik'
    options
    fitErr                % fitting error
    ordgpMdl              % an OrdRegressionGP object
    covFcnType = {'squaredexponential', ...
                  'ardsquaredexponential'} % covariance functions accepted by OrdRegressionGP
    logModel              % display model object after each training
    binEdges              % edges of binning
  end

  methods (Access = public)
    function obj = OrdGpModel(modelOptions, xMean)
      % constructor
      assert(size(xMean,1) == 1, 'OrdGpModel (constructor): xMean is not a row-vector.');

      % modelOpts structure
      if (isempty(modelOptions))
        obj.options = struct();
      else
        obj.options = modelOptions;
      end

      % computed settings
      obj.dim       = size(xMean, 2);
      obj.shiftMean = zeros(1, obj.dim);
      obj.shiftY    = 0;

      % fitting depends on fmincon
      if (~license('checkout', 'optimization_toolbox'))
        warning('OrdGpModel: Optimization Toolbox license not available. Model cannot be trained');
      end

      obj.covFcn = defopts(obj.options, 'covFcn', '{@covMaterniso, 5}');
      if isfield(obj.options, 'hyp')
        obj.hyp.ordreg = defopts(obj.options.hyp, 'ordreg', []);
        % suppose the input values are logarithmic
        obj.hyp.cov = defopts(obj.options.hyp, 'cov', log([0.5; 2]));
        obj.hyp.lik = defopts(obj.options.hyp, 'lik', log(0.01));
      else
        obj.hyp.ordreg = [];
        obj.hyp.cov = log([0.5; 2]);
        obj.hyp.lik = log(0.01);
      end
      
      if (exist(obj.covFcn, 'file') == 2)
      % string with name of an m-file function
        obj.covFcn  = str2func(covFcn);
      elseif (~any(ismember(obj.covFcn, obj.covFcnType)))
        % a function handle
        if (~isfield(obj.hyp, 'cov') || isempty(obj.hyp.cov))
          error('Hyperparameters must be specified for custom covariance functions');
        end
        obj.covFcn = eval(obj.covFcn);
      end
      
      % normalize options
      obj.options.normalizeY = defopts(obj.options, 'normalizeY', true);
      obj.options.normalizeX = defopts(obj.options, 'normalizeX', true);

      % prediction types:
      %
      % * avgord - average of ordinal responses weighted by predicted probabilities
      % * metric - metric predictions of the latent variable, i.e. w/o mapping into ordinal values
      obj.options.prediction = defopts(obj.options, 'prediction', 'avgord');
      
      % binning settings
      obj.options.binning = defopts(obj.options, 'binning', 'unipoints');
      if strcmp(obj.options.binning, 'none')
        % binning type 'none' does not need the number of bins
        obj.options.nBins = defopts(obj.options, 'nBins', 0);
      else
        obj.options.nBins = defopts(obj.options, 'nBins', 'mu + 1');
      end
      % normalize settings
      obj.options.normalizeY = defopts(obj.options, 'normalizeY', true);

      % general model prediction options
      obj.predictionType = defopts(modelOptions, 'predictionType', 'fValues');
      obj.transformCoordinates = defopts(modelOptions, 'transformCoordinates', true);

      % the rest of initial values
      obj.ordgpMdl = [];
      obj.trainLikelihood = Inf;
      obj.logModel = 1;
      obj.fitErr.mzoe = [];
      obj.fitErr.mae = [];
    end

    function nData = getNTrainData(obj)
      % returns the required number of data for training the model
      % TODO: *write this* properly according to dimension and
      %       covariance function set in options
      %       It should be nBins dependent. However, this value is at this
      %       moment unknown...
      nData = 3 * obj.dim;
    end

    function obj = trainModel(obj, X, y, xMean, generation)
      assert(size(xMean,1) == 1, '  OrdGpModel.train(): xMean is not a row-vector.');
      obj.trainMean = xMean;
      obj.dataset.X = X;
      obj.dataset.y = y;
      
      % check the training data
      if ( size(X, 1) < obj.getNTrainData ) || ( length(y) < obj.getNTrainData )
        fprintf(2, 'OrdGpModel.train(): Not enough data for model training.\n');
        obj.trainGeneration = -1;
        return
      elseif size(X, 1) ~= length(y)
        fprintf(2, 'OrdGpModel.train(): Different number of training points in X and y.\n');
        obj.trainGeneration = -1;
        return
      end
      
      mu = obj.stateVariables.mu;
      lambda = obj.stateVariables.lambda;
      nBins = myeval(obj.options.nBins);
      nTrain = size(obj.dataset.X, 1);

      % normalize y if specified or if large y-scale
      % (at least for binning)
      if (~obj.options.normalizeY ...
          && (max(y) - min(y)) > 1e4)
        fprintf(2, 'OrdGpModel.train(): Y-Normalization is switched ON for large Y-scale.\n');
        obj.options.normalizeY = true;
      end
      if (obj.options.normalizeY)
        obj.shiftY = mean(y);
        obj.stdY = std(y);
        yTrain = (y - obj.shiftY) / obj.stdY;
      else
        obj.shiftY = 0;
        obj.stdY = 1;
        yTrain = y;
      end
      
      % compute logarigthm of the input
%       if any(yTrain <= 0)
%         yTrain = yTrain - min(yTrain, [], 1) + eps;
%         obj.yTrainNeg = true;
%       else
%         obj.yTrainNeg = false;
%       end
%       yTrain = log(yTrain);
      
      % perform binning of fitness logarithm
      [yTrainBin, obj.binEdges] = binning(yTrain, nBins, obj.options.binning);

      % ordgp options
      ordgpOpts = { ...
        'FitMethod', 'exact', ...
        'Standardize', obj.options.normalizeX, ...
        'KernelFunction', obj.covFcn, ...
        'NumStartPoints', 2
      };

      if (isfield(obj.hyp, 'ordreg') && ~isempty(obj.hyp.ordreg))
        ordgpOpts(end+1:end+2) = {'PlsorParameters', obj.hyp.ordreg};
      end

      if (isfield(obj.hyp, 'cov') && ~isempty(obj.hyp.cov))
        ordgpOpts(end+1:end+2) = {'KernelParameters', obj.hyp.cov};
      end
      
      if (isfield(obj.options, 'covBounds') && ~isempty(obj.options.covBounds))
        ordgpOpts(end+1:end+2) = {'KernelBounds', obj.options.covBounds};
      end

      if (isfield(obj.hyp, 'lik') && ~isempty(obj.hyp.lik))
        ordgpOpts(end+1:end+2) = {'Sigma2', obj.hyp.lik};
      end
      
      if (isfield(obj.options, 'likBounds') && ~isempty(obj.options.likBounds))
        ordgpOpts(end+1:end+2) = {'Sigma2Bounds', obj.options.likBounds};
      end

      % tic;
      % train ordinal regression model
      obj.ordgpMdl = OrdRegressionGP(obj.dataset.X, yTrainBin, ordgpOpts);
      % fprintf('Toc: %.2f\n', toc);
      
      % check the model accuracy
      [yTrainPredict, ~, ~, ~, yTrainPredict_exp] = obj.ordgpMdl.predict(obj.dataset.X);
      % mean absolute error
      obj.fitErr.mae = sum(abs(yTrainPredict - yTrainBin)) / nTrain;
      % mean absolute error of probability-weighted prediction
      obj.fitErr.maew = sum(abs(yTrainPredict_exp - yTrainBin)) / nTrain;
      % mean zero-one error
      obj.fitErr.mzoe = sum(yTrainPredict ~= yTrainBin) / nTrain;
      % fprintf('MAE: %0.4f  MAEW: %0.4f  MZOE: %0.4f\n', obj.fitErr.mae, obj.fitErr.maew, obj.fitErr.mzoe)
      
      if (obj.ordgpMdl.MinimumNLP < Inf) % && (obj.fitErr.mae <= nBins / 5)
        obj.trainLikelihood = obj.ordgpMdl.MinimumNLP;
        obj.trainGeneration = generation;
      else
        obj.trainLikelihood = Inf;
        obj.trainGeneration = -1;
      end

%       if (obj.logModel)
%         disp('Model:');
%         disp(obj.ordgpMdl);
%       end
    end

    function [ypred, ysd2] = modelPredict(obj, X)
      % predicts the function values in new points X
      % @ypred      -- predicted response
      % @ysd2       -- predicted variance
      if (strcmpi(class(obj.ordgpMdl), 'OrdRegressionGP'))
        [~, ~, ymu, gp_ysd2, y_exp] = obj.ordgpMdl.predict(X);

        % number of outputs
        n = size(ymu, 1);
        % number of bins
        nBins = length(obj.binEdges) - 1;

        if (n <= 0)
          warning('Empty prediction.');
          ypred = [];
          ysd2 = [];
          return;
        end
        
        switch obj.options.prediction
          case 'avgord'
            ypred = y_exp;
          case 'metric'
            alpha = obj.ordgpMdl.PlsorParameters(1);
            ypred = (alpha < 0)*nBins + sign(alpha)*ymu;
          otherwise
            ypred = ymu;
        end
        
        % normalize minimal values
        ymin_norm = (min(obj.dataset.y) - obj.shiftY) / obj.stdY;
        ymax_norm = (max(obj.dataset.y) - obj.shiftY) / obj.stdY;
        % prediction correction using original binning edges
        newEdges = [ymin_norm, obj.binEdges(2:end-1), ymax_norm]';
        % find apropriate index in binning edges
        ymu_int = min(max(1, floor(ypred + 1/2)), length(newEdges) - 1);
        ymu_rem = ypred - ymu_int + 1/2;
        % scale the first bin according to minimal ymu
        ymu_rem(ypred < 1/2) = 1/2*(ymu_rem(ypred < 1/2) - min(ymu_rem))/(1/2 - min(ymu_rem));
        ypred = (newEdges(ymu_int + 1) - newEdges(ymu_int)).*ymu_rem + newEdges(ymu_int);
        
        % un-normalize in the f-space
        ypred = ypred * obj.stdY + obj.shiftY;
        ysd2 = gp_ysd2 * (obj.stdY)^2;
        
      else
        ypred = [];
        ysd2 = [];
        warning('Model not trained');
      end
    end

    function mdl = getModel(obj)
      % get OrdRegressionGP object
      mdl = obj.ordgpMdl;
    end
    
    function edges = getBinEdges(obj)
      % get binning edges of model
      edges = obj.binEdges;
    end
    
    function trainLik = getTrainLikelihood(obj)
      % get training likelihood
      trainLik = obj.trainLikelihood;
    end
    
    function opts = getOptions(obj)
      % get model options
      opts = obj.options;
    end
    
    function res = myeval(s)
      if ischar(s)
        res = evalin('caller', s);
      else
        res = s;
      end
    end

  end % methods

end % classdef
