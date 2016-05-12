function [xbest, fmin, counteval, stopflag, y_eval] = gpop(fitfun, xstart, gpopOpts, cmOpts, modelOpts)
% Implementation of Gaussian Process Optimization Procedure.
%
% Procedure uses a Gaussian process as a surrogate model
% of the objective function that is trained and optimized
% with CMA-ES algorithm at each iteration.
%
% The process is trained on most recently evaluated data
% from currently best point's neighbourhood.
%
% @fitfun                 -- objective function to minimize
% @xtart                  -- starting point
% @gpopOpts
%   'nc'                  -- number of training points selected by euclidean distance from 'xbest'
%   'nr'                  -- number of training points selected by evaluation time
%   'meritParams'         -- parameters that define surrogate 'merit' functions
%   'tolXPrtb'            -- tolerance for x improvement before a perturbation is added to 'xbest'
%   'prtb'                -- standard deviation of Gaussian perturbation added on stagnation
%   'maxFunEvals'         -- maximum number of fitness evaluations
%   'maxIter'             -- maximum number of iterations
%   'maxIterPrtb'         -- maximum number of consecutive stagnating iterations
%   'stopFitness'         -- stop if fitness exceeds specified value, minimization
%   'tolFunHist'          -- stop if range of recorded fitness changes smaller than 'tolHistFun'
%   'funHistLen'          -- length of recorded fitness history
%   'logModulo'           -- log inner variables every nth iteration
% @cmOpts                 -- CMA-ES options
% @modelOpts              -- GP model options
%
% Reference: D. Buche, N. N. Schraudolph and P. Koumoutsakos,
% "Accelerating Evolutionary Algorithms With Gaussian Process Fitness Function Models",
% IEEE Transactions on Systems, Man, and Cybernetics Part C (Applications and Reviews),
% 06/2005, 35(2):183-194.

dim = size(xstart, 1);

% default options
opts = struct(...
  'nc', '5 * dim', ...
  'nr', '5 * dim', ...
  'meritParams', [0 1 2 4], ...
  'tolXPrtb', 1e-8, ...
  'prtb', 1e-2, ...
  'maxFunEvals', 1e4, ...
  'maxIter', '1e3 * dim^2 / sqrt(opts.nc + opts.nr)', ...
  'maxIterPrtb', 'opts.nc + opts.nr', ...
  'stopFitness', -Inf, ...
  'tolFunHist', 1e-9, ...
  'funHistLen', 2, ...
  'logModulo', 100 ...
);

for fname = fieldnames(gpopOpts)'
  opts.(fname{1}) = gpopOpts.(fname{1});
end

xbest = xstart;
fmin = Inf;
fhist = NaN(opts.funHistLen, 1);
countiter = 0;
counteval = 0;
stopflag = {};
iterPrtb = 0;
y_eval = [];
surrogateStats = NaN(1, 2);

% eval string parameters
if ischar(modelOpts.hyp.cov), modelOpts.hyp.cov = eval(modelOpts.hyp.cov); end

for fname = fieldnames(opts)'
  if ischar(opts.(fname{1})), opts.(fname{1}) = eval(opts.(fname{1})); end
end

stopflag = stop_criteria();

if isempty(stopflag)
  countiter = countiter + 1;

  model = GpopGpModel(modelOpts, ones(1, dim));
  archive = GpopArchive(dim);

  % initial search for nc / 2 points with (2,10)-CMA-ES
  cmOptsInit = cmOpts;
  cmOptsInit.ParentNumber = 2;
  cmOptsInit.PopSize = 10;
  cmOptsInit.MaxFunEvals = ceil(opts.nc / 2);

  [xCm1, fminCm1, countevalCm1, stopflagCm1, outCm1, besteverCm1, y_evalCm1] = s_cmaes(fitfun, xstart, 8/3, cmOptsInit);

  xbest = besteverCm1.x;
  fmin = besteverCm1.f;
  counteval = counteval + countevalCm1;
  y_eval = [y_eval; y_evalCm1];
  archive = archive.save(outCm1.arxvalids', outCm1.fvalues', countiter);
  stopflag = stop_criteria();
  log_state();
end

while isempty(stopflag)
  % % uncomment when debugging
  % try
  countiter = countiter + 1;
  sol = []; % array with solutions found for model's current state

  % select training data
  [closestX, closestY] = archive.getNearData(opts.nc, xbest');
  [recentX, recentY] = archive.getRecentData(opts.nr);
  trainingX = union(closestX, recentX, 'rows');
  trainingY = union(closestY, recentY);

  % train model
  model = model.trainModel(trainingX, trainingY, xbest', countiter);

  if model.isTrained()
    % restrict CMA-ES search area to xbest's neighbourhood
    d = max(closestX)' - min(closestX)';
    cmOpts.LBounds = xbest - d/2;
    cmOpts.UBounds = xbest + d/2;
    sigma = [];

    % optimize all variants of model prediction with parallel workers
    parfor i = 1:length(opts.meritParams)
      % minimize surrogate function
      fun = @(x) surrogateFcn(x, opts.meritParams(i), model);
      [xCm, fminCm, ~, stopflagCm, ~, besteverCm, ~] = s_cmaes(fun, xbest, sigma, cmOpts);
      res(:, i) = besteverCm.x;
    end % parfor

    % evaluate model optima and save new solutions to archive
    for x = res
      if (~archive.isMember(x', opts.tolXPrtb))
        y = eval_fitness(x);
        sol = [sol [x; y]];
        stopflag = stop_criteria();
        if ~isempty(stopflag)
          y_eval = [y_eval; fmin counteval surrogateStats];
          log_state();
          return;
        end % if
      end % if
    end % for
  end % if

  if isempty(sol)
    % no solution found or model not trained
    % add a Gaussian perturbation of xbest to the training data
    d = max(closestX)' - min(closestX)';
    x = xbest + opts.prtb * (randn(size(xbest)) .* d);
    y = eval_fitness(x);
    fhist(2:end) = fhist(1:end-1);
    fhist(1) = y;
    iterPrtb = iterPrtb + 1;
  else
    % save only the best solution in current iteration to the history
    fhist(2:end) = fhist(1:end-1);
    fhist(1) = min(sol(2,:));
    iterPrtb = 0;
  end

  stopflag = stop_criteria();

  y_eval = [y_eval; fmin counteval surrogateStats];

  if ~isempty(stopflag) || (mod(countiter-1, opts.logModulo) == 0)
    log_state();
  end
  % catch err
  %   disp(err.message);
  % end % catch
end % while

  function flag = stop_criteria()
    % return flags of any stopping criteria that hold true
    flag = {};
    if iterPrtb >= opts.maxIterPrtb, flag(end+1) = {'maxiterprtb'}; end
    if counteval >= opts.maxFunEvals, flag(end+1) = {'maxfunevals'}; end
    if countiter >= opts.maxIter, flag(end+1) = {'maxiter'}; end
    if fmin <= opts.stopFitness, flag(end+1) = { 'fitness' }; end
    if (countiter >= max([length(fhist) opts.nc])) && ...
        (max(fhist) - min(fhist) <= opts.tolFunHist)
      flag(end+1) = { 'tolfunhist' };
    end
  end % function

  function y = eval_fitness(x)
    y = feval(fitfun, x);

    counteval = counteval + 1;
    archive = archive.save(x', y, countiter);

    if y < fmin
      xbest = x;
      fmin = y;
    end % if
  end % function

  function log_state()
    varnames = { 'countiter', 'fmin', 'xbest', 'counteval', 'fchange', ...
      'iterPrtb', 'stopflag' };
    if isempty(stopflag)
      flag = { [] };
    else
      flag = stopflag;
    end
    t = table([countiter], [fmin], { mat2str(xbest', 2) }, [counteval], ...
      [max(fhist) - min(fhist)], iterPrtb, flag, ...
      'VariableNames', varnames);
    disp(t);
  end % function
end % function


function y = surrogateFcn(x, a, model)
  [mu, s2] = model.modelPredict(x');
  y = fmerit(mu, s2, a);
end % function


function y = fmerit(mu, s2, a)
  % merit function that adds a density measure to the predicted value
  % @mu       -- column vector of predicted means
  % @s2       -- column vector of predicted output variances
  % @a        -- scale parameter
  y = mu - a * sqrt(s2);
end % function
