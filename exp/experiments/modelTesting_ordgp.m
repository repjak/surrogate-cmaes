
% ordgp model testing

% settings
func = (1:24);
dims = [2, 5, 10, 20];
maxEvals = 100;

% path settings
defFolder = fullfile('exp', 'experiments', 'model');
defSetFile = fullfile(defFolder, 'defData', ['defSet_', num2str(maxEvals), 'FE.mat']);
defModelFolder = fullfile(defFolder, 'defData', ['defModel_', num2str(maxEvals), 'FE']);

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
modelTestSets('exp_doubleEC_21_log', func, dims, maxEvals);
% load dataset
% ds_load = load(fullfile(defFolder, 'defSet'));
% ds = ds_load;

%% test chosen models
modelFolders = testModels('ordgp', modelOptions, defSetFile, func, dims);
% add default folder (GpModel)
modelFolders = [defModelFolder; modelFolders];

%% compare results
modelStatistics(modelFolders, func, dims)