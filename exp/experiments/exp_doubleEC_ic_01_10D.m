exp_id = 'exp_doubleEC_ic_01';
exp_description = 'Adaptive DTS, model selection with information criteria, 7 GP models, 5 insts, 10D';

% BBOB/COCO framework settings

bbobParams = { ...
  'dimensions',         { 10 }, ...
  'functions',          num2cell(1:24), ...      % all functions: num2cell(1:24)
  'opt_function',       { @opt_s_cmaes }, ...
  'instances',          num2cell(1:5), ...    % default is [1:5, 41:50]
  'maxfunevals',        { '250 * dim' }, ...
  'resume',             { false }, ...
};

% Surrogate manager parameters

surrogateParams = { ...
  'evoControl',         { 'doubletrained' }, ...    % 'none', 'individual', 'generation', 'restricted'
  'observers',          { {'DTScreenStatistics', 'DTFileStatistics'} },... % logging observers
  'modelType',          { 'bayesmodelsel' }, ...               % 'gp', 'rf', 'bbob'
  'updaterType',        { 'rankDiff' }, ...         % OrigRatioUpdater
  'DTAdaptive_updateRate',     { 0.3 }, ...
  'DTAdaptive_updateRateDown', { 'obj.updateRate' }, ...
  'DTAdaptive_maxRatio',       { 1.0 }, ...
  'DTAdaptive_minRatio',       { 0.04 }, ...
  'DTAdaptive_lowErr',         { '@(x) [ones(size(x,1),1) log(x(:,1)) x(:,2) log(x(:,1)).*x(:,2) x(:,2).^2] * [0.11; -0.0092; -0.13; 0.044; 0.14]' }, ...
  'DTAdaptive_highErr',        { '@(x) [ones(size(x,1),1) log(x(:,1)) x(:,2) log(x(:,1)).*x(:,2) x(:,2).^2] * [0.35; -0.047; 0.44; 0.044; -0.19]' }, ...
  'DTAdaptive_defaultErr',     { 0.05 }, ...
  'evoControlMaxDoubleTrainIterations', { 1 }, ...
  'evoControlPreSampleSize',       { 0.75 }, ...       % {0.25, 0.5, 0.75}, will be multip. by lambda
  'evoControlNBestPoints',         { 0 }, ...
  'evoControlValidationGenerationPeriod', { 4 }, ...
  'evoControlValidationPopSize',   { 0 }, ...
  'evoControlOrigPointsRoundFcn',  { 'ceil' }, ...  % 'ceil', 'getProbNumber'
  'evoControlIndividualExtension', { [] }, ...      % will be multip. by lambda
  'evoControlBestFromExtension',   { [] }, ...      % ratio of expanded popul.
  'evoControlTrainRange',          { 10 }, ...      % will be multip. by sigma
  'evoControlTrainNArchivePoints', { '15*dim' },... % will be myeval()'ed, 'nRequired', 'nEvaluated', 'lambda', 'dim' can be used
  'evoControlSampleRange',         { 1 }, ...       % will be multip. by sigma
  'evoControlOrigGenerations',     { [] }, ...
  'evoControlModelGenerations',    { [] }, ...
  'evoControlValidatePoints',      { [] }, ...
  'evoControlRestrictedParam',     { 0.05 }, ...
};

% Hyperparameter priors
pg1 = '{@priorGauss, log(0.01), 1}';
pt1 = '{@priorT, log(1), 4, 5}';

% Specification of model set for model selector

modelOptions = struct( ...
  'name',   {'SE', 'NN', 'LIN', 'QUAD', 'ADD', 'SE+NN', 'SE+LIN'}, ...
  'type',   repmat({'gp'}, 1, 7), ...
  'params', { ...
    struct( ...   % SE
      'covFcn',   '@covSEiso', ...
      'hyp',      struct('lik', log(0.01), 'cov', log([1 1])), ...
      'prior',    struct('lik', {{pg1}}, 'cov', {{pt1 pt1}}) ...
    ), ...
    struct( ...   % NN
      'covFcn',   '@covNNone', ...
      'hyp',      struct('lik', log(0.01), 'cov', log([1 1])), ...
      'prior',    struct('lik', {{pg1}}, 'cov', {{pt1 pt1}}) ...
    ), ...
    struct( ...   % LIN
      'covFcn',   '@covLINone', ...
      'hyp',      struct('lik', log(0.01), 'cov', log(1)), ...
      'prior',    struct('lik', {{pg1}}, 'cov', {{pt1}}) ...
    ), ...
    struct( ...   % QUAD
      'covFcn',   '{@covPoly, ''eye'', 2}', ...
      'hyp',      struct('lik', log(0.01), 'cov', log([1 1])), ...
      'prior',    struct('lik', {{pg1}}, 'cov', {{pt1 pt1}}) ...
    ), ...
    struct( ...   % ADD
      'covFcn',   '{@covADD, {[1], @covSEisoU}}', ...
      'hyp',      struct('lik', log(0.01), 'cov', log([ones(1, 10) 1])), ... % CAUTION: dim-dependent!
      'prior',    struct('lik', {{pg1}}, 'cov', {[repmat({pt1}, 1, 10) pt1]}) ...
    ), ...
    struct( ...   % SE + NN
      'covFcn',   '{@covSum, {@covSEiso, @covNNone}}', ...
      'hyp',      struct('lik', log(0.01), 'cov', log([1 1 1 1])), ...
      'prior',    struct('lik', {{pg1}}, 'cov', {{pt1 pt1 pt1 pt1}}) ...
    ), ...
    struct( ...   % SE + LIN
      'covFcn',   '{@covSum, {@covSEiso, @covLINone}}', ...
      'hyp',      struct('lik', log(0.01), 'cov', log([1 1 1])), ...
      'prior',    struct('lik', {{pg1}}, 'cov', {{pt1 pt1 pt1}}) ...
    ) ...
  } ...
);

% Model parameters

modelParams = { ...
  'modelOptions', { modelOptions }, ...
  'sharedModelOptions', { ...
    struct( ...
      'meanFcn',            { 'meanZero' }, ...
      'likFcn',             { 'likGauss' }, ...
      'predictionType',     { 'poi' }, ...
      'useShift',           { false }, ...
      'normalizeY',         { true }, ...
      'centerX',            { true }, ...
      'trainAlgorithm',     { 'mcmc' }, ...
      'mcmcNSimu',          { 'max(100, 2 * obj.nHyp ^ 2)' }, ...
      'mcmcBurnin',         { 'obj.mcmcNSimu' }, ...
      'mcmcNChains',        { 2 }, ...
      'predictFullPost',    { true }, ...
      'nSimuPost',          { '10 * obj.nHyp' } ...
    ) ...
  }, ...
  'ic',                 { 'waic2' }, ...
  'factory',            { 'ModelFactory' }, ...
  'trainsetType',       { 'nearest' }, ... % inherited properties
  'predictionType',     { 'poi' }, ...
  'trainRange',         { 4 }, ...
  'trainsetSizeMax',    { '20*dim' }, ...
  'transformCoordinates', { true } ...
};

% CMA-ES parameters

cmaesParams = { ...
  'PopSize',            { '(8 + floor(6*log(N)))' }, ...        %, '(8 + floor(6*log(N)))'};
  'Restarts',           { 50 }, ...
  'DispModulo',         { 0 }, ...
};

logDir = '/storage/plzen1/home/repjak/public';