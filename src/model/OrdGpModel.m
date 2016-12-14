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
    fitErr
    ordgpMdl              % an OrdRegressionGP object
    covFcnType = {'squaredexponential', ...
      'ardsquaredexponential'} % covariance functions accepted by OrdRegressionGP
    logModel              % display model object after each training
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
      obj.hyp.ordreg = defopts(obj.options.hyp, 'ordreg', []);
      obj.hyp.cov = defopts(obj.options.hyp, 'cov', []);
      obj.hyp.lik = defopts(obj.options.hyp, 'lik', []);

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

      obj.ordgpMdl = [];
      obj.logModel = 1;
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

      % normalize y if specified or if large y-scale
      % (at least for CMA-ES hyperparameter optimization)
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
      
      % perform binning
      yTrain = binning(yTrain, 5);

      % ordgp options
      ordgpOpts = { ...
        'FitMethod', 'exact', ...
        'Standardize', obj.options.normalizeX, ...
        'KernelFunction', obj.options.covFcn
      };

      if (isfield(obj.hyp, 'ordreg') && ~isempty(obj.hyp.ordreg))
        ordgpOpts(end+1:end+2) = {'PlsorParameters', obj.hyp.ordreg};
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
      
      obj.trainGeneration = generation;

      if (obj.logModel)
        disp('Model:');
        disp(obj.ordgpMdl);
      end
    end

    function [ypred, ysd2] = modelPredict(obj, X)
      % make prediction
      % @ypred      -- predicted response
      % @ysd        -- predicted standard deviation
      if (strcmpi(class(obj.ordgpMdl), 'OrdRegressionGP'))
        [~, yprob, ymu, ysd2] = obj.ordgpMdl.predict(X);
        % ysd = sqrt(ysd2);

        % un-normalize in the f-space
        ymu = ymu * obj.stdY + obj.shiftY;
        % TODO: correct ysd2
        % ysd = ysd * obj.stdY;

        % number of outputs
        n = size(yprob, 1);

        % number of ordinal classes
        r = size(yprob, 2);

        if (n <= 0)
          warning('Empty prediction.');
          ypred = [];
          ysd2 = [];
          return;
        end

        switch obj.options.prediction
          case 'avgord'
            ypred = sum(bsxfun(@times, yprob, 1:r), 2) ./ n;
          case 'metric'
            ypred = ymu;
          otherwise
            ypred = ymu;
        end
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

  end % methods

end % classdef
