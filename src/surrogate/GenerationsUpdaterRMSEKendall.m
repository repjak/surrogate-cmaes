classdef GenerationsUpdaterRMSEKendall < GenerationsUpdater
  properties
    parsedParams
    minModelGenerations
    maxModelGenerations
    errThreshold
    updateRate
    rmse
    kendall
    err
    generations
    lastModelGenerations
    origGenerations
    steepness
    alpha
    transferFun
  end

  methods (Static)
    function y = simplesig(x, k)
      % a simple "normalized" sigmoid function obtaining values 0 and 1 in
      % points 0 and 1, respectively and having the middle point in 0.5
      % @k -- a steepness factor, k > 0
      r = 1 + 1/k;
      y = r * (x - 0.5) ./ (r - 1 + abs(2 * (x - 0.5))) + 0.5;
    end

    function Kmn = covFcn(Xm, Xn, theta)
      % A covariance function defined as sum of SEard, constant
      % shift and a white noise:
      %
      % Kmn(x_i, x_j) = sf2 * exp(-(x_i - x_q)'*diag(1./ell)*(x_i - x_j)/2)
      %                 + sc2
      %                 + sn2*delta(x_i, x_j)
      %
      % where ell = (ell_1^2,...ell_D^2) is vector of ARD parameters, D is
      % dimension of input space, sf2 is signal variance, sc2 is constant shift
      % and sn2 is added white noise.
      %
      % Theta is a vector of positive numbers:
      % theta = [ log(ell_1)
      %           log(ell_2)
      %            .
      %           log(ell_D)
      %           log(sqrt(sf2))
      %           log(sqrt(sc2))
      %           log(sqrt(sn2)) ]
      %
      assert(size(Xm, 2) == size(Xn, 2), 'Dimensions must agree');
      dim = size(Xm, 2);

      % tolerance on squared distance of vectors when their considered equal
      tol = 1e-9;

      assert(all(size(theta) == [dim+3 1]), ['theta must be a row vector of ''' dim + 3 ''' hyperparameters']);

      ell = exp(theta(1:dim));
      sf2 = exp(2*theta(dim+1));
      sc2 = exp(2*theta(dim+2));
      sn2 = exp(2*theta(dim+3));

      Kmn = sq_dist(diag(1./ell) * Xm', diag(1./ell) * Xn');
      Kmn = sf2 * exp(-Kmn/2);
      Kmn = bsxfun(@plus, Kmn, sc2);

      delta = bsxfun(@lt, sq_dist(Xm', Xn'), tol^2);
      Kmn(delta) = Kmn(delta) + sn2;
    end
  end

  methods
    function [newOrigGenerations, newModelGenerations] = update(obj, modelY, origY, ~, ~, countiter)
      obj.generations(end+1) = countiter;

      if isempty(modelY)
        obj.rmse(end+1) = NaN;
        obj.kendall(end+1) = NaN;
        obj.err(end+1) = NaN;
      else
        newRMSE = sqrt(sum((modelY - origY).^2))/length(origY);
        newKendall = corr(modelY, origY, 'type', 'Kendall');

        if ~isempty(obj.rmse) && ~all(isnan(obj.rmse))
          maxRMSE = max(obj.rmse);
        else
          maxRMSE = newRMSE;
        end

        % combine the RMSE and Kendall ranking coefficient
        newErr = (1 - obj.alpha) * (newRMSE / maxRMSE) + obj.alpha * (1 - newKendall) / 2;

        % exponential smoothing
        if ~isempty(obj.err) && ~all(isnan(obj.err))
          errValid = obj.err(~isnan(obj.err));
          lastErr = errValid(end);
        else
          lastErr = 0.5;
        end

        newErr = (1 - obj.updateRate) * lastErr + obj.updateRate * newErr;

        obj.rmse(end+1) = newRMSE;
        obj.kendall(end+1) = newKendall;
        obj.err(end+1) = newErr;

        % normalize into an interval below an error threshold and apply a
        % transfer function
        errNormalized = obj.transferFun(max(1 - newErr / obj.errThreshold, 0));

        % min-max scaling and rounding
        obj.lastModelGenerations = round(obj.minModelGenerations + errNormalized * (obj.maxModelGenerations - obj.minModelGenerations));

        disp(['GenerationsUpdaterRMSE: model generations: ', num2str(round(obj.lastModelGenerations)), ' [ ', repmat('+', 1, round(obj.lastModelGenerations)), repmat(' ', 1, obj.maxModelGenerations - round(obj.lastModelGenerations)), ' ]']);
      end

      newOrigGenerations = round(obj.origGenerations);
      newModelGenerations = round(obj.lastModelGenerations);
    end

    function obj = GenerationsUpdaterRMSEKendall(parameters)
      % constructor
      obj = obj@GenerationsUpdater(parameters);
      if ~isstruct(parameters)
        obj.parsedParams = struct(parameters{:});
      else
        obj.parsedParams = parameters;
      end
      obj.minModelGenerations = defopts(obj.parsedParams, 'minModelGenerations', 1);
      obj.maxModelGenerations = defopts(obj.parsedParams, 'maxModelGenerations', 5);
      obj.updateRate = defopts(obj.parsedParams, 'updateRate', 0.5);
      obj.errThreshold = defopts(obj.parsedParams, 'errThreshold', 0.45);
      obj.steepness = defopts(obj.parsedParams, 'steepness', 5);
      obj.alpha = defopts(obj.parsedParams, 'alpha', 0.2);
      obj.transferFun = defopts(obj.parsedParams, 'transferFun', @(x) x);
      obj.rmse = [];
      obj.err = [];
      obj.origGenerations = 1;
      obj.lastModelGenerations = obj.minModelGenerations;
    end
  end
end