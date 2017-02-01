
% ordgp model testing

% settings
func = (1:24);
dims = [2, 5, 10, 20];
maxEvals = 100;

% path settings
defFolder = fullfile('exp', 'experiments', 'model');
defModelFolder = fullfile(defFolder, 'defData', 'defModel');

% default model options
defModelOptions.useShift = false;
defModelOptions.predictionType = 'sd2';
defModelOptions.trainAlgorithm = 'fmincon';
defModelOptions.covFcn = '{@covMaterniso, 5}';
defModelOptions.normalizeY = true;
defModelOptions.hyp.lik = log(0.01);
defModelOptions.hyp.cov = log([0.5; 2]);

% test model options
modelSet = defModelOptions;
% none binning
modelSet.prediction = {'avgord', 'metric'};
modelSet.binning = 'none';
modelOptions = combineFieldValues(modelSet);
% the rest of settings
modelSet.binning = {'logcluster', 'unipoints'};
modelSet.nBins = {'mu', 'lambda', '2*lambda'};
modelRest = combineFieldValues(modelSet);
% add rest to model options
modelOptions = [modelOptions; modelRest];

%% create testing dataset
ds = modelTestSets('exp_doubleEC_21_log', func, dims, maxEvals);

%% test chosen models
modelFolders = testModels('ordgp', modelOptions, ds, func, dims);
% add default folder (GpModel)
modelFolders = [modelFolders; defModelFolder];

%% compare results
compareModels(modelFolders, func, dims)