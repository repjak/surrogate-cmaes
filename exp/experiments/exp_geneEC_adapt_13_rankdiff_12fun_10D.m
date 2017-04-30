exp_id = 'exp_geneEC_adapt_13_rankdiff_12fun_10D';
exp_description = 'Surrogate CMA-ES model using generation EC with generation number dependant on model''s RankDiff error on 12 functions, ardmatern52, gpml';

% BBOB/COCO framework settings

bbobParams = { ...
  'dimensions',       { 2, 5, 10 }, ...
  'functions',        { 2, 3, 6, 8, 12, 13, 15, 17, 18, 21, 23, 24 }, ... %num2cell(1:24), ...
  'opt_function',     { @opt_s_cmaes }, ...
  'instances',        { [1:5]; [41:45]; [41:50] }, ...
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
  'updaterType',                   { 'RankDiff' }, ...     % { 'rmseKendall', 'constant', 'KullbackLeibler' }
  'geneECAdaptive_minModelGenerations', { 0 }, ...
  'geneECAdaptive_maxModelGenerations', { 5, 10, 20 }, ...
  'geneECAdaptive_updateRate',     { 0.2, 0.5, 0.8 }, ...
  'geneECAdaptive_transferFun',    { '@(x)x', '@(x)GenerationsUpdater.simplesig(x, 2)' }, ...
  'geneECAdaptive_aggregateType',  { 'lastvalid' }, ...
  'geneECADaptive_lowErrThreshold', {0}, ...
  'geneECAdaptive_highErrThreshold', {0.5, 0.9} ...
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
