function [x, fmin, counteval, stopflag] = gpop(fitfun, xstart, gpopOpts, cmOpts, modelOpts)
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
%   'prtb'                -- standard deviation of Gaussian perturbation
%   'stopMaxIter'         -- maximum number of iterations
%   'stopMaxIterPrtb'     -- maximum number of consecutive iterations with little improvement
%   'stopMaxFunEvals'     -- maximum number of evaluations of @fitfun
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
  'tolXPrtb', 1e-3, ...
  'prtb', 1e-2, ...
  'stopMaxFunEvals', Inf, ...
  'stopMaxIterPrtb', 10, ...
  'stopMaxIter', '1e3 * dim^2 / sqrt(opts.nc + opts.nr)');

for fname = fieldnames(gpopOpts)'
  opts.(fname{1}) = gpopOpts.(fname{1});
end

x = Inf(2, 1);
xbest = xstart';
fmin = Inf;
countiter = 1;
counteval = 0;
stopflag = {};
iterPrtb = 0;

% eval string parameters
if ischar(modelOpts.hyp.cov)
  modelOpts.hyp.cov = eval(modelOpts.hyp.cov);
end

for fname = fieldnames(opts)'
  if ischar(opts.(fname{1}))
    opts.(fname{1}) = eval(opts.(fname{1}));
  end
end

model = GpModel(modelOpts, ones(1, dim));
archive = GpopArchive(dim);

% initial search for nc / 2 points with (2,10)-CMA-ES
cmOptsInit = cmOpts;
cmOptsInit.ParentNumber = 2;
cmOptsInit.PopSize = 10;
cmOptsInit.MaxFunEvals = ceil(opts.nc / 2);

[x1 fmin1 counteval1 stopflag1 out1 bestever1 y_eval1] = s_cmaes(fitfun, xstart, 8/3, cmOptsInit);

xbest = bestever1.x';
fmin = bestever1.f;

archive = archive.save(out1.arxvalids', out1.fvalues', countiter);

while isempty(stopflag)
  [closestX, closestY, closestIdx] = archive.getNearData(opts.nc, xbest);
  % diameter of the smallest hypercube containing all closest point
  diam = max(closestX) - min(closestX);
  [recentX, recentY, recentIdx] = archive.getRecentData(opts.nr, xbest, diam, closestIdx);
  trainingX = [closestX; recentX];
  trainingY = [closestY; recentY];

  % alternatively: don't filter previously selected points in the second step
  % this may result in training sets smaller than nc + nr
  % [recent, recentY] = archive.getRecentData(opts.nr, xbest, diam);
  % trainingX = union(closestX, recentX, 'rows');
  % trainingY = union(closestY, recentY, 'rows');

  % train model
  model = model.trainModel(trainingX, trainingY, xbest, countiter);

  % optimize all variants of model prediction
  for a = opts.meritParams
    sol = {}; % array with solutions found in current iteration

    % minimize surrogate function
    fun = @(x) surrogateFcn(x, a, model);
    [xCm fminCm countevalCm stopflagCm outCm besteverCm y_evalCm] = s_cmaes(fun, xbest, 8/3, cmOpts);
    counteval = counteval + countevalCm;

    % save new solutions to archive
    if (~archive.isMember(besteverCm.x', opts.tolXPrtb))
      sol(end+1) = {besteverCm};
      archive = archive.save(besteverCm.x', besteverCm.f, countiter);
    end
  end % for

  if isempty(sol)
    % add a Gaussian perturbation to the training data
    iterPrtb = iterPrtb + 1;
    x = xbest + randn() * diam * opts.prtb;
    y = fun(x');
    counteval = counteval + 1;
    archive = archive.save(x, y, countiter);
  else
    % update xbest
    for i = 1:length(sol)
      if sol{i}.f < fmin
        xbest = sol{i}.x';
        fmin = sol{i}.f;
      end
    end
    iterPrtb = 0;
  end

  % stop criteria
  if iterPrtb >= opts.stopMaxIterPrtb, stopflag(end+1) = {'maxiterprtb'}; end
  if counteval >= opts.stopMaxFunEvals, stopflag(end+1) = {'maxfunevals'}; end
  if countiter >= opts.stopMaxIter, stopflag(end+1) = {'maxiter'}; end

  countiter = countiter + 1;
  x = xbest';
end % while

end % function


function y = surrogateFcn(x, a, model)
  [mu, s2] = model.modelPredict(x');
  if (model.options.normalizeY)
    % un-normalize predicted deviation
    s2 = s2 / model.stdY ^ 0.5;
  end
  y = fmerit(mu, s2, a);
end % function


function y = fmerit(mu, s2, a)
  % merit function that adds a density measure to the predicted value
  % @mu       -- column vector of predicted means
  % @s2       -- column vector of predicted output variances
  % @a        -- scale parameter
  y = mu - a * sqrt(s2);
end % function
