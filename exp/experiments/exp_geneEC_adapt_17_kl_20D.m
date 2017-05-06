exp_id = 'exp_geneEC_adapt_17_kl_20D';
exp_description = 'Surrogate CMA-ES model using generation EC with generation number dependant on model''s KullbackLeibler divergence, best settings, ardmatern52';

% BBOB/COCO framework settings

bbobParams = { ...
  'dimensions',       { 20 }, ...
  'functions',        num2cell(1:24), ...
  'opt_function',     { @opt_s_cmaes }, ...
  'instances',        num2cell([1:5, 41:50]), ...
  'maxfunevals',      { '250 * dim' }, ...
};

% Surrogate manager parameters

surrogateParams = { ...
  'evoControl',                    { 'generation' }, ...
  'modelType',                     { 'gp' }, ...       % 'gp', 'fitrgp', 'rf', 'bbob'
  'evoControlPreSampleSize',       { 0 }, ...              % {0.25, 0.5, 0.75}, will be multip. by lambda
  'evoControlTrainRange',          { 8 }, ...              % will be multip. by sigma
  'evoControlTrainNArchivePoints', { '15*dim' }, ...       % will be myeval()'ed, 'nRequired', 'nEvaluated', 'lambda', 'dim' can be used
  'evoControlSampleRange',         { 1 }, ...              % will be multip. by sigma
  'evoControlValidatePoints',      { [] }, ...
  'updaterType',                   { 'KullbackLeibler' }, ...     % { 'rmseKendall', 'constant', 'RankDiff' }
  'geneECAdaptive_minModelGenerations', { 0 }, ...
  'geneECAdaptive_maxModelGenerations', { 5 }, ...
  'geneECAdaptive_updateRate',     { 0.5 }, ...
  'geneECAdaptive_transferFun',    { '@(x)GenerationsUpdater.simplesig(x, 2)' }, ...
  'geneECAdaptive_aggregateType',  { 'lastvalid' }, ...
  'geneECADaptive_lowErrThreshold', { 0 }, ...
  'geneECAdaptive_highErrThreshold', { 0.9 } ...
};

% Model parameters

modelParams = { ...
  'useShift',           { false }, ...
  'trainAlgorithm',     { 'fmincon' }, ...
  'covFcn',             { '{@covMaterniso, 5}' }, ...
  'hyp',                { struct('lik', log(0.01), 'cov', log([0.5; 2])) }, ...
  'normalizeX',         { true }, ...
  'normalizeY',         { true } ...
};

% CMA-ES parameters

cmaesParams = { ...
  'PopSize',            { '(4 + floor(3*log(N)))' }, ...
  'Restarts',           { 4 }, ...
};

logDir = '/storage/plzen1/home/repjak/public';
