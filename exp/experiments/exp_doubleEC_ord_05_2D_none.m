exp_id = 'exp_doubleEC_ord_05_2D_none';
exp_description = 'Surrogate CMA-ES model, ordinal GP, covSE, no clustering, doubly-trained EC, 24 functions, 15 instances, 2D, (source code revised)';

% BBOB/COCO framework settings

bbobParams = { ...
  'dimensions',         { 2 }, ...
  'functions',          num2cell(1:24), ...      % all functions: num2cell(1:24)
  'opt_function',       { @opt_s_cmaes }, ...
  'instances',          { 1:5, 41:45, 46:50 }, ...    % default is [1:5, 41:50]
  'maxfunevals',        { '100 * dim' }, ...
  'resume',             { true }
};

% Surrogate manager parameters

surrogateParams = { ...
  'evoControl',         { 'doublytrained' }, ...    % 'none', 'individual', 'generation', 'restricted'
  'observers',          { {'DTScreenStatistics', 'DTFileStatistics'} },... % logging observers
  'modelType',          { 'ordgp' }, ...            % 'gp', 'ordgp', 'rf', 'bbob'
  'evoControlPreSampleSize', { 0 }, ...             % {0.25, 0.5, 0.75}, will be multip. by lambda
  'evoControlIndividualExtension', { [] }, ...      % will be multip. by lambda
  'evoControlBestFromExtension', { [] }, ...        % ratio of expanded popul.
  'evoControlTrainRange', { 10 }, ...               % will be multip. by sigma
  'evoControlTrainNArchivePoints', { '15*dim' },... % will be myeval()'ed, 'nRequired', 'nEvaluated', 'lambda', 'dim' can be used
  'evoControlSampleRange', { 1 }, ...               % will be multip. by sigma
  'evoControlOrigGenerations', { [] }, ...
  'evoControlModelGenerations', { [] }, ...
  'evoControlValidatePoints', { [] }, ...
  'evoControlRestrictedParam', { 0.05 }, ...
};

% Model parameters

modelParams = { ...
  'useShift',           { false }, ...
  'trainAlgorithm',     { 'fmincon' }, ...
  'covFcn',             { 'squaredexponential' }, ...
  'nBestPoints',        { 0 }, ...
  'minLeaf',            { 2 }, ...
  'inputFraction',      { 1 }, ...
  'normalizeY',         { true }, ...
  'prediction',         { 'avgord'}, ...
  'binning',            { 'none' }, ...
  'nBins',              { 0 }
};

% CMA-ES parameters

cmaesParams = { ...
  'PopSize',            { '(8 + floor(6*log(N)))' }, ...        %, '(8 + floor(6*log(N)))'};
  'Restarts',           { 4 }, ...
};

logDir = '/storage/plzen1/home/bajeluk/public';
