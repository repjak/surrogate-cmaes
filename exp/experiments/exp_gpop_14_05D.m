exp_id = 'exp_gpop_14_05D';
exp_description = 'GPOP (Buche) with GP and CMA-ES, 3 funs, 3 instances, 100 * dim maxfunevals, 3 various sets of sigma and hyperparameters for custom cov fcn, slightly lower funHistLen';

% BBOB/COCO framework settings

bbobParams = { ...
  'dimensions',              { 5 }, ...
  'functions',               { 1 }, ...
  'opt_function',            { @opt_gpop }, ...
  'instances',               { [1 2 41] }, ...     % default is [1:5, 41:50]
  'maxfunevals',             { '100 * dim' }, ...
};

% GPOP parameters

surrogateParams = { ...
  'gpop_nc',                 { '5 * dim' }, ...       % number of training points selected by distance to xbest
  'gpop_nr',                 { 'opts.nc' }, ...       % number of training points selected by time of evaluation
  'gpop_meritParams',        { [0 1 2 4] }, ...       % parameters of merit function (added density measure)
  'gpop_tolXPrtb',           { 1e-8 }, ...            % tolerance for improvement before a perturbation is considered
  'gpop_maxIterPrtb',        { 2 }, ...               % maximum number of consecutive perturbations
  'gpop_tolFunHist',         { 1e-9 }, ...            % stop if range of recorded f-values lower than tolFunHist
  'gpop_funHistLen',         { 5 }, ...               % length of recorded f-values history
};

% Model parameters

modelParams = { ...
  'cov',                     { struct( ...                                                       % settings we used with gpml
                                 'fcn', '@covFcn', ...
                                 'sigma', 0.001, ...
                                 'hyp', 'log([0.5 * ones(dim, 1); 0.1; 0.1; 1e-4])' ...
                               ), ...
                               struct( ...
                                 'fcn', '@covFcn', ...
                                 'sigma', 0.001, ...
                                 'hyp', 'log([5 * ones(dim, 1); 0.5; 0.5; 0.05])' ...
                               ), ...
                               struct( ...
                                 'fcn', '@covFcn', ...
                                 'sigma', 0.001, ...
                                 'hyp', 'log([5 * ones(dim, 1); 0.5; 0.5; 1e-4])' ...
                               ), ...
                             }, ...
  'normalizeY',              { true }, ...
  'normalizeX',              { true } ...
};

%                               struct( ...
%                                 'fcn', '@covFcn', ...
%                                 'sigma', 0.01, ...
%                                 'hyp', 'log([5 * ones(dim, 1); 0.5; 0.5; 0.05])' ...
%                               ), ...

% CMA-ES parameters

cmaesParams = { ...
  'PopSize',                 { 10 }, ...
  'ParentNumber',            { 2 }, ...
  'Restarts',                { 4 }, ...
};

logDir = '/storage/plzen1/home/repjak/public';
