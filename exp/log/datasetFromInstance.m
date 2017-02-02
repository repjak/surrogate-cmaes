function dataset = datasetFromInstance(exp_id, nDatasetsPerInstance, fun, dim, id, maxFunEvals)
%DATASETFROMINSTANCE - generates datasets for offline model tunning
%
% Generates datasets for offline model tunning from the file with 
% saved models (*_modellog_*.mat) and DTS-CMA-ES results file
% (*_results_*.mat)

  % Initialization

  % DEBUG:
  % exp_id = 'exp_doubleEC_21_log';
  % nDatasetsPerInstance = 10;
  % fun = 1;
  % dim = 2;
  % id = 1;

  dataset = {};

  if nargin < 6
    if nargin < 1
      help datasetFromInstance
      return
    end
    maxFunEvals = 250;
  end
  
  experimentPath = [pwd '/exp/experiments/' exp_id];

  % load data from files
  savedModelsFile = sprintf('%s/bbob_output/%s_modellog_%d_%dD_%d.mat', experimentPath, exp_id, fun, dim, id);
  scmaesOutFile = sprintf('%s/%s_results_%d_%dD_%d.mat', experimentPath, exp_id, fun, dim, id);
  if exist(savedModelsFile, 'file') && exist(scmaesOutFile, 'file')
    MF = load(savedModelsFile, 'models');
    models = MF.models;
    SF = load(scmaesOutFile, 'cmaes_out', 'exp_settings');
    cmaes_out = SF.cmaes_out;
    exp_settings = SF.exp_settings;
  else
    warning('Model file or scmaes output file is missing in f%d %dD (id %d).', fun, dim, id)
    return
  end

  % BBOB fitness initialization
  fgeneric('initialize', exp_settings.bbob_function, exp_settings.instances(1), ['/tmp/bbob_output/']);

  % find maximal evaluation generation
  maxGener = find(cmaes_out{1}{1}.origEvaled, maxFunEvals*dim, 'first');
  nGenerations = cmaes_out{1}{1}.generations(maxGener(end));
  gens = floor(linspace(0, nGenerations-1, nDatasetsPerInstance+1));
  gens(1) = [];
  cmo = cmaes_out{1}{1};

  trainSetX = cell(nDatasetsPerInstance, 1);
  trainSetY = cell(nDatasetsPerInstance, 1);
  testSetX = cell(nDatasetsPerInstance, 1);
  testSetY = cell(nDatasetsPerInstance, 1);

  % return no result if there is not enough models
  % TODO: solve this
  if length(models) < gens(end)
    warning('Not enough models to generate dataset!')
    return
  end
  
  % Dataset generation

  for i = 1:nDatasetsPerInstance
    g = gens(i);
    m = models{g};

    lambda = sum(cmo.generations == g);

    fprintf('Train set size (gen.# %3d): %3d\n', g, size(m.dataset.X,1));

    cmaesState = struct( ...
      'xmean', m.trainMean', ...
      'sigma', m.trainSigma, ...
      'lambda', lambda, ...
      'BD', m.trainBD, ...
      ... % 'diagD', diagD, ...
      'diagD', [], ...
      ... % 'diagC', diagC, ...
      'dim', m.dim, ...
      'mu', floor(lambda/2), ...
      'countiter', g);

    sampleOpts = struct( ...
      'noiseReevals', 0, ...
      'isBoundActive', true, ...
      'lbounds', -5 * ones(m.dim, 1), ...
      'ubounds',  5 * ones(m.dim, 1), ...
      'counteval', cmaes_out{1}{1}.generationStarts(g), ...
      'flgEvalParallel', false, ...
      'flgDiagonalOnly', false, ...
      'noiseHandling', false, ...
      'xintobounds', @xintobounds, ...
      'origPopSize', lambda);

    % Generate fresh CMA-ES' \lambda offsprings
    [arx, arxvalid, arz] = sampleCmaesNoFitness(m.trainSigma, lambda, cmaesState, sampleOpts);

    % Save everything needed
    trainSetX{i} = (m.trainBD * m.dataset.X' * m.trainSigma)';
    trainSetY{i} = m.dataset.y;
    testSetX{i}  = arxvalid';
    testSetY{i}  = fgeneric(testSetX{i}')';
    means{i}     = m.trainMean;
    sigmas{i}    = m.trainSigma;
    BDs{i}       = m.trainBD;
    cmaesStates{i} = cmaesState;
  end

  % Finalize

  dataset = struct();
  dataset.trainSetX = trainSetX;
  dataset.trainSetY = trainSetY;
  dataset.testSetX  = testSetX;
  dataset.testSetY  = testSetY;
  dataset.means = means;
  dataset.sigmas = sigmas;
  dataset.BDs = BDs;
  dataset.cmaesStates = cmaesStates;
  dataset.generations = gens;
  dataset.function  = fun;
  dataset.dim       = m.dim;
  dataset.id        = id;
  dataset.sampleOpts = sampleOpts;

  fgeneric('finalize');
end