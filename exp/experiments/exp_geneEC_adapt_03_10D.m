exp_id = 'exp_geneEC_adapt_03_10D';
exp_description = 'Surrogate CMA-ES model using generation EC with generation number dependant on model''s log-scaled RMSE and Kendall correlation';

% BBOB/COCO framework settings

bbobParams = { ...
  'dimensions',        { 2, 3, 5, 10 }, ...
  'functions',        num2cell(1:24), ...
  'opt_function',     { @opt_s_cmaes }, ...
  'instances',        { [1:5, 41:50] }, ...
  'maxfunevals',      { '250 * dim' }, ...
};

% Surrogate manager parameters

surrogateParams = { ...
  'evoControl',                    { 'generation' }, ...
  'modelType',                     { 'gp' }, ...           % 'gp', 'rf', 'bbob'
  'evoControlPreSampleSize',       { 0 }, ...              % {0.25, 0.5, 0.75}, will be multip. by lambda
  'evoControlIndividualExtension', { 100 }, ...            % will be multip. by lambda
  'evoControlBestFromExtension',   { 0.2 }, ...            % ratio of expanded popul.
  'evoControlTrainRange',          { 8 }, ...              % will be multip. by sigma
  'evoControlTrainNArchivePoints', { '15*dim' }, ...       % will be myeval()'ed, 'nRequired', 'nEvaluated', 'lambda', 'dim' can be used
  'evoControlSampleRange',         { 1 }, ...              % will be multip. by sigma
  'evoControlValidatePoints',      { [] }, ...
  'updaterType',                   { 'rmse' }, ...         % { 'rmse', 'constant' }
  'updaterParams',                 { struct( ...
                                       'updaterType', 'rmse', ...
                                       'steepness',  5, ...
                                       'minModelGenerations', 1, ...
                                       'maxModelGenerations', 10 ...
                                     ), ...
                                     struct( ...
                                       'updaterType', 'rmse', ...
                                       'steepness',  5, ...
                                       'minModelGenerations', 1, ...
                                       'maxModelGenerations', 5 ...
                                     ), ...
                                   }, ...
};

% Model parameters

modelParams = { ...
  'useShift',           { false }, ...
  'trainAlgorithm',     { 'fmincon' }, ...
  'covFcn',             { '{@covMaterniso, 5}' }, ...
  'hyp',                { struct('lik', log(0.01), 'cov', log([0.5; 2])) }, ...
  'nBestPoints',        { 0 }, ...
  'minLeaf',            { 2 }, ...
  'inputFraction',      { 1 }, ...
  'normalizeY',         { true }, ...
};

% CMA-ES parameters

cmaesParams = { ...
  'PopSize',            { '(4 + floor(3*log(N)))' }, ...
  'Restarts',           { 4 }, ...
};

logDir = '/storage/plzen1/home/repjak/public';
