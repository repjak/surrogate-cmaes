
% ordgp model testing

% settings
func = 1;
dims = 2;

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

% create testing dataset
ds = modelTestSets('exp_doubleEC_21_log', func, dims);

% test chosen models
modelFolders = testModels('ordgp', modelOptions, ds, func, dims);

% compare results
compareModels(modelFolders, func, dims)