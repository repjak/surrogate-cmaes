classdef GprModel < Model
  % Implements GP with RegressionGP from MATLAB
  properties    % inherited from abstract class "Model"
    dim                   % dimension of the input space X (determined from x_mean)
    trainGeneration = -1; % # of the generation when the model was built
    trainMean             % mean of the generation when the model was trained
    trainSigma            % sigma of the generation when the model was trained
    trainBD               % BD of the generation when the model was trained
    dataset               % .X and .y
    useShift = false;
    shiftMean             % vector of the shift in the X-space
    shiftY = 0;           % shift in the f-space
    predictionType        % type of prediction (f-values, PoI, EI)
    transformCoordinates  % transform X-space
    stateVariables        % variables needed for sampling new points as CMA-ES do
    
    % GprModel specific properties
    stdY                  % standard deviation of Y in training set, for output normalization
    covFcn                % covariance function passed to fitrgp
    options
    hyp
    nErrors
    gprMdl                % a RegressionGP object
  end

  properties (Access = protected)
    % a list of strings identifiers of covariance functions supported by
    % fitrgp
    covFcnType = {'squaredexponential', 'matern32', 'matern52', 'ardsquaredexponential', ...
      'ardmatern32', 'ardmatern52'}
  end

  methods (Access = public)
    function obj = GprModel(modelOptions, xMean)
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

      % Statistics and ML Toolbox check
      if (~license('checkout', 'statistics_toolbox'))
        warning('GprModel: Statistics and Machine Learning Toolbox license not available. Model cannot be used');
      end

      % Optimization Toolbox check
      obj.options.trainAlgorithm = defopts(obj.options, 'trainAlgorithm', 'fmincon');
      if (strcmpi(obj.options.trainAlgorithm, 'fmincon') ...
          && ~license('checkout', 'optimization_toolbox'))
        warning('GpModel: Optimization Toolbox license not available. Switching to minimize().');
        obj.options.trainAlgorithm = 'fminsearch';
      end

      covFcn = defopts(obj.options, 'covFcn',  '{@covMaterniso, 5}');
      if (any(ismember(covFcn, obj.covFcnType)))
        % a string
        obj.covFcn = covFcn;
      else
        % a function handle
        obj.covFcn = eval(covFcn);
      end

      dim = obj.dim;

      obj.options.normalizeY = defopts(obj.options, 'normalizeY', true);
      obj.options.normalizeX = defopts(obj.options, 'normalizeX', true);
      obj.transformCoordinates = defopts(modelOptions, 'transformCoordinates', true);
      obj.hyp.sigma = defopts(obj.options.hyp, 'sigma', 'std(y)/sqrt(2)');
      obj.hyp.cov = defopts(obj.options.hyp, 'cov', '[std(X), std(y)/sqrt(2)]''');
      disp(obj.hyp.cov);
      obj.gprMdl = [];
    end

    function nData = getNTrainData(obj)
      nData = 3 * obj.dim;
    end
    
    function trained = isTrained(obj)
      trained = (strcmpi(class(obj.gprMdl), 'RegressionGP'));
    end

    function obj = trainModel(obj, X, y, xMean, generation)
      obj.dataset.X = X;
      obj.dataset.y = y;

      % normalize y if specified
      if (obj.options.normalizeY)
        obj.shiftY = mean(y);
        obj.stdY = std(y);
        y = (y - obj.shiftY) / obj.stdY;
      else
        obj.shiftY = 0;
        obj.stdY = 1;
      end

      alg = obj.options.trainAlgorithm;
      hyp = myeval(obj.hyp.cov);
      sigma = myeval(obj.hyp.sigma);

%       if (strcmpi(alg, 'fmincon') && (~isfield(obj.options, 'lb_hyp') ...
%           || ~isfield(obj.options, 'ub_hp')))
%         warning('Traning algorithm ''fmincon'' specified but no constraints given. Falling back to MATLAB''s default');
%         alg = 'default';
%       end

      if (strcmpi(alg, 'fminsearch'))
        obj.gprMdl = fitrgp(obj.dataset.X, y, ...
          'FitMethod', 'exact', ...
          'Sigma', sigma, ...
          'KernelFunction', obj.covFcn, ...
          'KernelParameters', hyp, ...
          'Optimizer', obj.options.trainAlgorithm, ...
          'Standardize', obj.options.normalizeX ...
        );
      elseif (strcmpi(alg, 'fmincon'))
        % TODO: find a way to pass lb, ub
        obj.gprMdl = fitrgp(obj.dataset.X, y, ...
          'FitMethod', 'exact', ...
          'PredictMethod', 'exact', ...
          'Sigma', sigma, ...
          'KernelFunction', obj.covFcn, ...
          'KernelParameters', hyp, ...
          'Optimizer', obj.options.trainAlgorithm, ...
          'Standardize', obj.options.normalizeX ...
        );
      else
        error(['Training algorithm ''' alg ''' not supported']);
      end

      disp('Model:');
      disp(obj.gprMdl);

    end

    function [ypred, ysd] = modelPredict(obj, X)
      % make prediction
      % @ypred      -- predicted response
      % @ysd        -- predicted standard deviation
      if (strcmpi(class(obj.gprMdl), 'RegressionGP'))
        [ypred, ysd] = predict(obj.gprMdl, X);
        ypred = ypred * obj.stdY + obj.shiftY;
        ysd = ysd * obj.stdY;
      else
        ypred = [];
        ysd = [];
        warning('Model not trained');
      end
    end
    
    function gprMdl = getGpr(obj)
      % get gprMdl object
      gprMdl = obj.gprMdl;
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
