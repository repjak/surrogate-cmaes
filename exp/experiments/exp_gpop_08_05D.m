exp_id = 'exp_gpop_08_05D';
exp_description = 'GPOP (Buche) with GP and CMA-ES, Sphere, Rastrigin, Rosenbrock, 2 instances, 100 * dim maxfunevals -- one task per instance';

% BBOB/COCO framework settings

bbobParams = { ...
  'dimensions',              { 5 }, ...
  'functions',               num2cell([2, 3, 8]), ...       % 2 Sphere, 3 Rastrigin, 8 Rosenbrock
  'opt_function',            { @opt_gpop }, ...
  'instances',               num2cell([1 41]), ...     % default is [1:5, 41:50]
  'maxfunevals',             { '100 * dim' }, ...
};

% GPOP parameters

surrogateParams = { ...
  'gpop_nc',                 { '5 * dim' }, ...      % number of training points selected by distance to xbest
  'gpop_nr',                 { 'opts.nc' }, ...      % number of training points selected by time of evaluation
  'gpop_meritParams',        { [0 1 2 4] }, ...      % parameters of merit function (added density measure)
  'gpop_tolXPrtb',           { 1e-8 }, ...           % tolerance for improvement before a perturbation is considered
  'gpop_maxIterPrtb',        { 2 }, ...              % maximum number of consecutive perturbations
  'gpop_tolFunHist',         { 1e-9 }, ...           % stop if range of recorded f-values lower than tolFunHist
  'gpop_funHistLen',         { 5 }, ...              % length of recorded f-values history
  'gpop_prtb',               { 1e-2 }, ...
  'gpop_logModulo',          { 1 }, ...
};

% Model parameters

modelParams = { ...
  'trainAlgorithm',          { 'gpop' }, ...
  'covFcn',                  { '{@covSum, {@covSEard, @covConst, @covNoise}}' }, ...
  'hyp',                     { struct( ...
                                 'lik', log(0.001), ...
                                 'cov', 'log([0.5 * ones(dim, 1); 0.1; 0.01; 1e-4])' ... % ell, signal variance, const shift, signal noise
                               ), ...
                             }, ...
  'normalizeY',              { true }, ...
};

% CMA-ES parameters

cmaesParams = { ...
  'PopSize',                 { 10 }, ...
  'ParentNumber',            { 2 }, ...
  'Restarts',                { 4 }, ...
};

logDir = '/storage/plzen1/home/repjak/public';
