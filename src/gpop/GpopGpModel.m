classdef GpopGpModel < GpModel
  % A GP model with training more in according with the article.

  properties
    nFminconIt;
    maxFminconIt;
    shiftX;
    stdX;
  end

  methods (Access = public)
    function obj = GpopGpModel(modelOptions, xMean)
      obj@GpModel(modelOptions, xMean);
      obj.nFminconIt = -1;  % track number of consecutive iterations with fmincon training
      obj.maxFminconIt = 20;   % max number of consecutive fmincon iterations before CMA-ES is used and counter reset
    end

    function obj = trainModel(obj, X, y, xMean, generation)
      if (~strcmpi(obj.options.trainAlgorithm, 'gpop'))
        obj = trainModel@GpModel(obj, X, y, xMean, generation);
        return
      end

      assert(size(xMean,1) == 1, '  GpModel.train(): xMean is not a row-vector.');
      obj.trainMean = xMean;
      obj.dataset.X = X;
      obj.dataset.y = y;

      % normalize y if specified, @meanZero, or if large y-scale
      % (at least for CMA-ES hyperparameter optimization)
      if (~obj.options.normalizeY ...
          && (isequal(obj.meanFcn, @meanZero) || (max(y) - min(y)) > 1e4))
        fprintf(2, 'Y-Normalization is switched ON for @meanZero covariance function of large Y-scale.\n');
        obj.options.normalizeY = true;
      end
      if (obj.options.normalizeY)
        obj.shiftY = mean(y);
        obj.stdY  = std(y);
        yTrain = (y - obj.shiftY) / obj.stdY;
      else
        obj.shiftY = 0;
        obj.stdY  = 1;
        yTrain = y;
      end

      % normalize X
      obj.shiftX = mean(X);
      obj.stdX = std(X);
      xTrain = (X - repmat(obj.shiftX, size(X, 1), 1)) ./ repmat(obj.stdX, size(X, 1), 1);

      % set the mean hyperparameter if is needed
      if (~isequal(obj.meanFcn, @meanZero))
        obj.hyp.mean = median(yTrain);
      end

      % gp() with linearized version of the hyper-parameters
      f = @(par) linear_gp(par, obj.hyp, obj.infFcn, obj.meanFcn, obj.covFcn, obj.likFcn, obj.dataset.X, yTrain);

      linear_hyp = unwrap(obj.hyp)';
      l_cov = length(obj.hyp.cov);

      assert(l_cov == 8, 'GPOP assumes exactly 8 hyperparameters');

      % lower and upper bounds
      lb_hyp.cov = log([1e-2 * ones(l_cov-3, 1); 1e-3; 1e-3; 1e-9]);
      ub_hyp.cov = log([10 * ones(l_cov-3, 1); 1; 1; 1e-2]);
      lb_hyp.lik = log(1e-6);
      ub_hyp.lik = log(10);
      % set bounds for mean hyperparameter
      if (~isequal(obj.meanFcn, @meanZero))
        minY = min(yTrain);
        maxY = max(yTrain);
        lb_hyp.mean = minY - 2*(maxY - minY);
        ub_hyp.mean = minY + 2*(maxY - minY);
      end

      lb = unwrap(lb_hyp)';
      ub = unwrap(ub_hyp)';

      opt = [];
      trainErr = false;

      if (obj.nFminconIt >= 0 && obj.nFminconIt < obj.maxFminconIt)
        [obj, opt, trainErr] = obj.trainFmincon(linear_hyp, xTrain, yTrain, lb, ub, f);

        if (trainErr)
          disp('Fmincon train error');
        else
          obj.nFminconIt = obj.nFminconIt + 1;
        end
      end

      if (obj.nFminconIt == -1 || obj.nFminconIt >= obj.maxFminconIt || trainErr)
        if (trainErr)
          disp('Trying CMA-ES...');
        end
        [obj, opt, trainErr] = obj.trainCmaes(linear_hyp, xTrain, yTrain, lb, ub, f);
        if (trainErr)
          % DEBUG OUTPUT:
          fprintf('.. model is not successfully trained, likelihood = %f\n', obj.trainLikelihood);
          return;
        end

        obj.nFminconIt = 0;
      end

      obj.trainGeneration = generation;
      obj.hyp = rewrap(obj.hyp, opt);

      % DEBUG OUTPUT:
      fprintf('.. model-training likelihood = %f\n', obj.trainLikelihood);
      % disp(obj.hyp);
    end

    function [y, dev] = modelPredict(obj, X)
      % predicts the function values in new points X
      if (obj.isTrained())
        % prepare the training set (if was normalized for training)
        nRows = size(obj.dataset.X, 1);
        xTrain = (obj.dataset.X - repmat(obj.shiftX, nRows, 1)) ./ repmat(obj.stdX, nRows, 1);
        yTrain = (obj.dataset.y - obj.shiftY) / obj.stdY;

        XNormalized = (X - repmat(obj.shiftX, size(X, 1), 1)) ./ repmat(obj.stdX, size(X, 1), 1);

        % apply the shift if the model is already shifted
        % TODO: do this before or after normalization in X?
        XWithShift = XNormalized - repmat(obj.shiftMean, size(XNormalized,1), 1);

        % calculate GP models' prediction in X
        [y, dev] = gp(obj.hyp, obj.infFcn, obj.meanFcn, obj.covFcn, obj.likFcn, xTrain, yTrain, XWithShift);
        % un-normalize in the f-space (if there is any)
        y = y * obj.stdY + obj.shiftY;
        dev = dev * obj.stdY;
      else
        y = []; dev = [];
        fprintf(2, 'GpModel.predict(): the model is not yet trained!\n');
      end
    end
  end

  methods (Access = private)
    function [obj, opt, trainErr] = trainFmincon(obj, linear_hyp, X, y, lb, ub, f);
      % train with Matlab's fmincon() from the Optimization toolbox
      %
      global modelTrainNErrors;
      trainErr = false;
      opt = [];

      [fminconOpts, nonlnc] = obj.defaultFminconOpts();
      try
        initial = f(linear_hyp');
      catch err
        initial = NaN;
      end
      if isnan(initial)
        % the initial point is not valid
        disp('  GpModel.train(): fmincon -- initial point is not valid.');
        trainErr = true;
      else
        % training itself
        disp(['Model training (fmincon), init fval = ' num2str(initial)]);
        try
          modelTrainNErrors = 0;
          [opt, fval] = fmincon(f, linear_hyp', [], [], [], [], lb, ub, nonlnc, fminconOpts);
          obj.nErrors = modelTrainNErrors;
          obj.trainLikelihood = fval;
          if (isnan(fval)  ||  initial - fval < 0.1)
            % final likelihood is not a valid value or
            % the shift in likelihood is almost none, the model is probably
            % not trained, do not use it
            trainErr = true;
          end
        catch err
          obj.nErrors = modelTrainNErrors;
          fprintf(2, '  GpModel.train() ERROR: fmincon() ended with an exception: %s\n', err.message);
          trainErr = true;
        end
      end
    end

    function [opts, nonlnc] = defaultFminconOpts(obj)
      % return the optimization parameters for fmincon()
      %
      opts = optimset('fmincon');
      opts = optimset(opts, ...
        'GradObj', 'on', ...
        'TolFun', 1e-6, ...
        'TolX', 1e-7, ...
        'MaxIter', 500, ...
        'MaxFunEvals', 500, ...
        'Display', 'off' ...
        );
      covarianceDim = length(obj.hyp.cov) - 1;
      if (covarianceDim > 1)
        % ARD
        opts = optimset(opts, 'Algorithm', 'interior-point');
        nonlnc = @nonlincons;
      else
        % ISOtropic
        opts = optimset(opts, 'Algorithm', 'trust-region-reflective');
        nonlnc = [];
      end
    end

    function [obj, opt, trainErr] = trainCmaes(obj, linear_hyp, X, y, lb, ub, f)
      % train with CMA-ES
      %
      global modelTrainNErrors;

      opt = []; fval = Inf; trainErr = false;
      cmaesopt.LBounds = lb';
      cmaesopt.UBounds = ub';
      cmaesopt.SaveVariables = false;
      cmaesopt.LogModulo = 0;
      cmaesopt.DispModulo = 0;
      cmaesopt.DispFinal = 0;
      sigma = [0.2*(ub - lb)]' + eps;
      % sigma(end) = min(10*mean(sigma(1:end-1)), sigma(end));
      if (length(obj.hyp.cov) > 2)
        % there is ARD covariance
        % try run cmaes for 500 funevals to get bounds for covariances
        MAX_DIFF = 2.5;
        cmaesopt.MaxFunEvals = 500;
        modelTrainNErrors = 0;
        try
          [opt, fval] = s_cmaes(f, linear_hyp', sigma, cmaesopt);
        catch err
          fprintf(2, 'GpModel.train() ERROR: CMA-ES ended with an exception: %s\n', err.message);
          trainErr = true;
          obj.nErrors = modelTrainNErrors;
          obj.trainGeneration = -1;
          return;
        end
        cov_median = median(opt(1:obj.dim));
        ub(1:obj.dim) = cov_median + MAX_DIFF;
        lb(1:obj.dim) = cov_median - MAX_DIFF;
        cmaesopt.LBounds = lb';
        cmaesopt.UBounds = ub';
        sigma(1:obj.dim) = [0.2*(ub(1:obj.dim) - lb(1:obj.dim))]' + eps;
      end
      cmaesopt.MaxFunEvals = 3000;
      try
        modelTrainNErrors = 0;
        [opt, fval] = s_cmaes(f, linear_hyp', sigma, cmaesopt);
      catch err
        fprintf(2, 'GpModel.train() ERROR: CMA-ES ended with an exception: %s\n', err.message);
        trainErr = true;
        obj.nErrors = modelTrainNErrors;
        obj.trainGeneration = -1;
        return;
      end
      if (isnan(fval))
        % final likelihood is not a valid value, the model is probably
        % not trained, do not use it
        trainErr = true;
      end
      obj.nErrors = modelTrainNErrors;
      obj.trainLikelihood = fval;
    end
  end % methods
end

function [nlZ, dnlZ] = linear_gp(linear_hyp, s_hyp, inf, mean, cov, lik, x, y)
  hyp = rewrap(s_hyp, linear_hyp');
  [nlZ, s_dnlZ] = gp(hyp, inf, mean, cov, lik, x, y);
  dnlZ = unwrap(s_dnlZ)';
end

function [c, ceq] = nonlincons(x)
  % checks if the values x(1:(end-4)) are within 2.5 off median
  MAX_DIFF = 2.5;
  ceq = [];
  assert(size(x,2) == 1, 'Argument for nonlincons is not a vector');
  c = zeros(size(x));
  % test only for covariance parameters
  % TODO: are there always 4 more parameters?!
  c = abs(x(1:end-4) - median(x(1:end-4))) - MAX_DIFF;
end
