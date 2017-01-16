classdef OrdGpModel < Model
  properties    % derived from abstract class "Model"
    dim                  % dimension of the input space X (determined from x_mean)
    trainGeneration      % # of the generation when the model was built
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

      obj.covFcn = defopts(obj.options, 'covfcn', 'ardsquaredexponential');
      if isfield(obj.options, 'hyp')
        obj.hyp.ordreg = defopts(obj.options.hyp, 'ordreg', []);
        obj.hyp.cov = defopts(obj.options.hyp, 'cov', []);
        obj.hyp.lik = defopts(obj.options.hyp, 'lik', []);
      else
        obj.hyp.ordreg = [];
        obj.hyp.cov = [];
        obj.hyp.lik = [];
      end

      if (~any(ismember(obj.covFcn, obj.covFcnType)))
        % a function handle
        if (~isfield(obj.hyp, 'cov') || isempty(obj.hyp.cov))
          error('Hyperparameters must be specified for custom covariance functions');
        end
        obj.covFcn = eval(obj.covFcn);
      end

      obj.options.normalizeY = defopts(obj.options, 'normalizeY', true);
      obj.options.normalizeX = defopts(obj.options, 'normalizeX', true);

      % prediction types:
      %
      % * avgord - average of ordinal responses weighted by predicted probabilities
      % * metric - metric predictions of the latent variable, i.e. w/o mapping into ordinal values
      obj.options.prediction = defopts(obj.options, 'prediction', 'avgord');
      obj.options.binning = defopts(obj.options, 'binning', 'unipoints');
      obj.options.nBins = defopts(obj.options, 'nBins', 'mu + 1');

      obj.ordgpMdl = [];
      obj.logModel = 1;
      obj.fitErr.mzoe = [];
      obj.fitErr.mae = [];
    end

    function nData = getNTrainData(obj)
      % returns the required number of data for training the model
      % TODO: *write this* properly according to dimension and
      %       covariance function set in options
      nData = 3 * obj.dim;
    end

    function obj = trainModel(obj, X, y, xMean, generation)
      assert(size(xMean,1) == 1, '  OrdGpModel.train(): xMean is not a row-vector.');
      obj.trainMean = xMean;
      obj.dataset.X = X;
      obj.dataset.y = y;
      
      mu = obj.stateVariables.mu;
      lambda = obj.stateVariables.lambda;
      nBins = myeval(obj.options.nBins);
      nTrain = size(obj.dataset.X, 1);

      % normalize y if specified or if large y-scale
      % (at least for binning)
      if (~obj.options.normalizeY ...
          && (max(y) - min(y)) > 1e4)
        fprintf(2, 'Y-Normalization is switched ON for large Y-scale.\n');
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
      
      % perform binning of fitness logarithm
      [yTrain, obj.binEdges] = binning(yTrain, nBins, obj.options.binning);

      % ordgp options
      ordgpOpts = { ...
        'FitMethod', 'exact', ...
        'Standardize', obj.options.normalizeX, ...
        'KernelFunction', obj.options.covFcn
      };

      if (isfield(obj.hyp, 'ordreg') && ~isempty(obj.hyp.ordreg))
        ordgpOpts(end+1:end+2) = {'PlsorParameters', obj.hyp.ordreg};
      % else 
      %  alpha = 1;
      %  beta1 = obj.binEdges(2);
      %  delta = diff(obj.binEdges(2:end-1));
      %  ordgpOpts(end+1:end+2) = {'PlsorParameters', [alpha, beta1, delta]};
      end

      if (isfield(obj.hyp, 'cov') && ~isempty(obj.hyp.cov))
        ordgpOpts(end+1:end+2) = {'KernelParameters', obj.hyp.cov};
      end

      if (isfield(obj.hyp, 'lik') && ~isempty(obj.hyp.lik))
        ordgpOpts(end+1:end+2) = {'Sigma2', obj.hyp.lik};
      end

      tic;
      % train ordinal regression model
      obj.ordgpMdl = OrdRegressionGP(obj.dataset.X, yTrain, ordgpOpts);
      fprintf('Toc: %.2f\n', toc);
      
      [yTrainPredict, ~, ~, ~, yTrainPredict_exp] = obj.ordgpMdl.predict(obj.dataset.X);
      % mean absolute error
      obj.fitErr.mae = sum(abs(yTrainPredict_exp - yTrain)) / nTrain;
      % mean zero-one error
      obj.fitErr.mzoe = sum(yTrainPredict ~= yTrain) / nTrain;
      fprintf('MAE: %0.4f  MZOE: %0.4f\n', obj.fitErr.mae, obj.fitErr.mzoe)
      
      if (obj.ordgpMdl.MinimumNLP < Inf) && (obj.fitErr.mae <= nBins / 5)
        obj.trainGeneration = generation;
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
        [~, yprob, ymu, gp_ysd2, y_exp] = obj.ordgpMdl.predict(X);

        % un-normalize in the f-space
%         ymu = ymu * obj.stdY + obj.shiftY;
%         ysd2 = gp_ysd2 * (obj.stdY)^2;
        ysd2 = gp_ysd2;

        % number of outputs
        n = size(yprob, 1);

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
            ypred = ymu;
          otherwise
            ypred = ymu;
        end
        
        % prediction correction using original binning edges
        newEdges = [min(obj.dataset.y), exp(obj.binEdges(2:end-1)), max(obj.dataset.y)]';
        % find apropriate index in binning edges
        ymu_int = min(max(1, floor(ypred + 1/2)), length(newEdges) - 1);
        ymu_rem = ypred - ymu_int + 1/2;
        % scale the first bin according to minimal ymu
        ymu_rem(ypred < 1/2) = 1/2*(ymu_rem(ypred < 1/2) - min(ymu_rem))/(1/2 - min(ymu_rem));
        ypred = (newEdges(ymu_int + 1) - newEdges(ymu_int)).*ymu_rem + newEdges(ymu_int);
        
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
    
    function res = myeval(s)
      if ischar(s)
        res = evalin('caller', s);
      else
        res = s;
      end
    end

  end % methods

end % classdef
