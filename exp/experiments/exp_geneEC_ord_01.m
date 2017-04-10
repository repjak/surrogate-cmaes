exp_id = 'exp_geneEC_ord_01';
exp_description = 'Surrogate CMA-ES model, ordinal GP, generation EC, 24 functions, 15 instances';

% BBOB/COCO framework settings

bbobParams = { ...
  'dimensions',         { 5 }, ...
  'functions',          { 2, 3, 6, 8, 13, 17, 21}, ...      % all functions: num2cell(1:24)
  'opt_function',       { @opt_s_cmaes }, ...
  'instances',          { [1:5, 41:50] }, ...    % default is [1:5, 41:50]
  'maxfunevals',        { '250 * dim' }, ...
  'resume',             { true }
};

% Surrogate manager parameters

surrogateParams = { ...
  'evoControl',         { 'generation' }, ...    % 'none', 'individual', 'generation', 'restricted'
  'observers',          { {'DTScreenStatistics', 'DTFileStatistics', 'DTModelSaver'} },... % logging observers
  'modelType',          { 'ordgp' }, ...               % 'gp', 'rf', 'bbob'
  'evoControlPreSampleSize', { 0 }, ...             % {0.25, 0.5, 0.75}, will be multip. by lambda
  'evoControlIndividualExtension', { [] }, ...      % will be multip. by lambda
  'evoControlBestFromExtension', { [] }, ...        % ratio of expanded popul.
  'evoControlTrainRange', { 10 }, ...               % will be multip. by sigma
  'evoControlTrainNArchivePoints', { '15*dim' },... % will be myeval()'ed, 'nRequired', 'nEvaluated', 'lambda', 'dim' can be used
  'evoControlSampleRange', { 1 }, ...               % will be multip. by sigma
  'evoControlOrigGenerations', { 1 }, ...
  'evoControlModelGenerations', { 1, 5 }
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
  'binning',            { 'none', 'logcluster', 'unipoints' }, ...
  'nBins',              { 'mu', 'lambda', '2*lambda' }
};

% CMA-ES parameters

cmaesParams = { ...
  'PopSize',            { '(4 + floor(3*log(N)))' }, ...        %, '(8 + floor(6*log(N)))'};
  'Restarts',           { 4 }, ...
};

logDir = '/storage/plzen1/home/bajeluk/public';