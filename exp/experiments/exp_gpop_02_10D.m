exp_id = 'exp_gpop_02_10D';
exp_description = 'GPOP (Buche) with GP and CMA-ES, 4 functions, all instances, max 10 * dim funevals';

% BBOB/COCO framework settings

bbobParams = { ...
  'dimensions',              { 2, 5, 10 }, ...
  'functions',               { 2, 3, 8, 20 } ...      % Buche, 2005: 2 Sphere, 3 Rastrigin, 8 Rosenbrock, 20 Schwefel; all: num2cell(1:24)
  'opt_function',            { @opt_gpop }, ...
  'instances',               { [1:5 41:50] }, ...     % default is [1:5, 41:50]
  'maxfunevals',             { '10 * dim' }, ...
  'parpoolSize',             { 4 }, ...
};

% GPOP parameters

surrogateParams = { ...
  'gpop_nc',                 { '5 * dim' }, ...       % number of training points selected by distance to xbest
  'gpop_nr',                 { 'opts.nc' }, ...       % number of training points selected by time of evaluation
  'gpop_meritParams',        { [0 1 2 4] }, ...       % parameters of merit function (added density measure)
  'gpop_tolXPrtb',           { 1e-8 }, ...            % tolerance for improvement before a perturbation is considered
  'gpop_maxIterPrtb',        { 2 }, ...               % maximum number of consecutive perturbations
  'gpop_tolFunHist',         { 1e-9 }, ...            % stop if range of recorded f-values lower than tolFunHist
  'gpop_funHistLen',         { 2 }, ...               % length of recorded f-values history
  };

% Model parameters

modelParams = { ...
  'trainAlgorithm',          { 'fmincon' }, ...
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
  'PopSize',                 { '(4 + floor(3*log(N)))', 10 }, ...
  'ParentNumber',            { 2 }, ...
  'Restarts',                { 4 }, ...
};

logDir = '/storage/plzen1/home/repjak/public';
