classdef ordgpMdlModel < Model
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
    function obj = ordgpMdlModel(modelOptions, xMean)
      % constructor
      assert(size(xMean,1) == 1, 'GpModel (constructor): xMean is not a row-vector.');

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
        warning('ordgpMdlModel: Optimization Toolbox license not available. Model cannot be trained');
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
      obj.logModel = 0;
    end

    function nData = getNTrainData(obj)
      % returns the required number of data for training the model
      % TODO: *write this* properly according to dimension and
      %       covariance function set in options
      nData = 3 * obj.dim;
    end

    function obj = trainModel(obj, X, y, xMean, generation)
      assert(size(xMean,1) == 1, '  GpModel.train(): xMean is not a row-vector.');

      obj.dataset.X = X;
      obj.dataset.y = y;

      % normalize y if specified or if large y-scale
      % (at least for CMA-ES hyperparameter optimization)
      if (~obj.options.normalizeY ...
          && (max(y) - min(y)) > 1e4))
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

      ordgpOpts = { ...
        'FitMethod', 'exact', ...
        'PredictMethod', 'exact', ...
        'Standardize', obj.options.normalizeX, ...
        'KernelFunction', obj.options.covFcn
      };

      if (isfield(obj.hyp.ordreg) && ~isempty(obj.hyp.ordreg))
        ordgOpts(end+1:end+2) = {'PlsorParameters', obj.hyp.ordreg};
      end

      if (isfield(obj.hyp.cov) && ~isempty(obj.hyp.cov))
        ordgpOpts(end+1:end+2) = {'KernelParameters', obj.hyp.cov};
      end

      if (isfield(obj.hyp.lik) && ~isempty(obj.hyp.lik))
        ordgpOpts(end+1:end+2) = {'Sigma2', obj.hyp.lik};
      end

      obj.ordgpMdl = OrdRegressionGP(obj.dataset.X, y, ordgpOpts);

      if (obj.logModel)
        disp('Model:');
        disp(obj.ordgpMdl);
      end
    end

    function [ypred, ysd] = modelPredict(obj, X)
      % make prediction
      % @ypred      -- predicted response
      % @ysd        -- predicted standard deviation
      if (strcmpi(class(obj.ordgpMdl), 'OrdRegressionGP'))
        [~, yprob, ymu, ysd2] = ordgpMdl.predict(X);
        ysd = sqrt(ysd2);

        % un-normalize in the f-space
        ypred = ypred * obj.stdY + obj.shiftY;
        ysd = ysd * obj.stdY;

        % number of outputs
        n = size(yprob, 1);

        % number of ordinal classes
        r = size(yprob, 2);

        if (n <= 0)
          warning('Empty prediction.');
          ypred = [];
          ysd = [];
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
        ysd = [];
        warning('Model not trained');
      end
    end

    function mdl = getModel(obj)
      % get OrdRegressionGP object
      mdl = obj.ordgpMdl;
    end

  end % methods

end % classdef
