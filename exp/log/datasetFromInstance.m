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
  models  = {};

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
  if exist(scmaesOutFile, 'file')
    if exist(savedModelsFile, 'file')
      MF = load(savedModelsFile, 'models');
      models = MF.models;
    end
    SF = load(scmaesOutFile, 'cmaes_out', 'exp_settings', 'surrogateParams');
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
  
  % Dataset generation

  for i = 1:nDatasetsPerInstance
    g = gens(i);

    lambda = sum(cmo.generations == g);

    if (length(models) >= g  &&  ~isempty(models{g}))
      % use the data in the model from DTS-CMA-ES run
      m = models{g};
      xmean = m.trainMean';
      sigma = m.trainSigma;
      BD    = m.trainBD;
      dim   = m.dim;
      X_train = (BD * m.dataset.X' * sigma)';
      y_train = m.dataset.y;
    else
      % the model is not saved, so create a fresh new one
      xmean = cmo.means(:,g);
      sigma = cmo.sigmas(g);
      BD    = cmo.BDs{g};
      dim   = exp_settings.dim;
      m = ModelFactory.createModel(SF.surrogateParams.modelType, SF.surrogateParams.modelOpts, xmean');

      % create the Archive of original points
      minTrainSize = m.getNTrainData();
      archive = Archive(dim);
      orig_id = logical(cmo.origEvaled(1:(cmo.generationStarts(g)-1)));
      X_train = cmo.arxvalids(:,orig_id)';
      y_train = cmo.fvalues(orig_id)';
      archive.save(X_train, y_train, 5);
      archive.gens = cmo.generations(orig_id);

      % find the model's trainset -- points near the current xmean
      nArchivePoints = myeval(SF.surrogateParams.evoControlTrainNArchivePoints);
      [X_train, y_train, nData] = archive.getDataNearPoint(nArchivePoints, ...
        xmean', SF.surrogateParams.evoControlTrainRange, ...
        sigma, BD);
    end

    fprintf('Train set size (gen.# %3d): %3d\n', g, size(X_train,1));

    cmaesState = struct( ...
      'xmean', xmean, ...
      'sigma', sigma, ...
      'lambda', lambda, ...
      'BD', BD, ...
      ... % 'diagD', diagD, ...
      'diagD', [], ...
      ... % 'diagC', diagC, ...
      'dim', dim, ...
      'mu', floor(lambda/2), ...
      'countiter', g);

    sampleOpts = struct( ...
      'noiseReevals', 0, ...
      'isBoundActive', true, ...
      'lbounds', -5 * ones(dim, 1), ...
      'ubounds',  5 * ones(dim, 1), ...
      'counteval', cmo.generationStarts(g), ...
      'flgEvalParallel', false, ...
      'flgDiagonalOnly', false, ...
      'noiseHandling', false, ...
      'xintobounds', @xintobounds, ...
      'origPopSize', lambda);

    % Generate fresh CMA-ES' \lambda offsprings
    [arx, arxvalid, arz] = sampleCmaesNoFitness(sigma, lambda, cmaesState, sampleOpts);

    % Save everything needed
    trainSetX{i} = X_train;
    trainSetY{i} = y_train;
    testSetX{i}  = arxvalid';
    testSetY{i}  = fgeneric(testSetX{i}')';
    means{i}     = xmean';
    sigmas{i}    = sigma;
    BDs{i}       = BD;
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