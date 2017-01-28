function overallLik = compareModels(modelOptions, dataset, comp_model, outputFolder, saveIter)
% overallLik = compareModels(modelOptions, dataset, comp_model, 
% outputFolder, saveIter) tests ordGpModel on datasets and compares results
% to GpModel.
%
% Input:
%   modelOptions - model options to test
%   dataset      - data for testing | Mx1 cell-array of string
%   comp_model   - model for testing according to dataset | Mx1 cell-array 
%                  of string
%   outputFolder - folder for saving results
%   saveIter     - save individual iterations (if compareModels is being
%                  optimized using external function) | boolean
%
% Output:
%   overallLik - sum of computed likelihoods on all datasets
%
% See Also:
%   OrdGpModel, GpModel

  defDataset = fullfile('exp', 'experiments', 'exp_doubleEC_21_short_log');
  if nargout > 0
    overallLik = Inf;
  end

  if nargin < 5
    if nargin < 4
      if nargin < 3
        if nargin < 2
          if nargin < 1
            help compareModels
            return
          end
          dataset = {fullfile(defDataset, 'exp_doubleEC_21_short_log_results_2_5D_1.mat'), ...
                     fullfile(defDataset, 'exp_doubleEC_21_short_log_results_3_5D_2.mat'), ...
                     fullfile(defDataset, 'exp_doubleEC_21_short_log_results_6_5D_3.mat'), ...
                     fullfile(defDataset, 'exp_doubleEC_21_short_log_results_8_5D_4.mat'), ...
                     fullfile(defDataset, 'exp_doubleEC_21_short_log_results_13_5D_5.mat'), ...
                     fullfile(defDataset, 'exp_doubleEC_21_short_log_results_17_5D_6.mat'), ...
                     fullfile(defDataset, 'exp_doubleEC_21_short_log_results_21_5D_7.mat')};
        end
        comp_model = {fullfile(defDataset, 'bbob_output', 'exp_doubleEC_21_short_log_modellog_2_5D_1.mat'), ...
                      fullfile(defDataset, 'bbob_output', 'exp_doubleEC_21_short_log_modellog_3_5D_2.mat'), ...
                      fullfile(defDataset, 'bbob_output', 'exp_doubleEC_21_short_log_modellog_6_5D_3.mat'), ...
                      fullfile(defDataset, 'bbob_output', 'exp_doubleEC_21_short_log_modellog_8_5D_4.mat'), ...
                      fullfile(defDataset, 'bbob_output', 'exp_doubleEC_21_short_log_modellog_13_5D_5.mat'), ...
                      fullfile(defDataset, 'bbob_output', 'exp_doubleEC_21_short_log_modellog_17_5D_6.mat'), ...
                      fullfile(defDataset, 'bbob_output', 'exp_doubleEC_21_short_log_modellog_21_5D_7.mat')};
      end
      outputFolder = 'exp/experiments/modelComparison';
    end
    saveIter = false;
  end
  
  % defaults
  nGen = 5;
  if ~isdir(outputFolder)
    mkdir(outputFolder)
  end
  
  % model options parsing
  if isnumeric(modelOptions) && ~isempty(modelOptions)
    fprintf('modelPoint: %s\n', num2str(modelOptions, '%0.4f  '))
    modelOpts = modelOptsParse(modelOptions);
  else
    modelOpts = modelOptions;
  end
  
  data = [];
  comparisonModel = [];
  
  % data input
  if isstruct(dataset)
    data = dataset;
  elseif ~iscell(dataset)
    dataset = {dataset};
  end
  nData = length(dataset);
  if iscell(comp_model) && iscell(comp_model{1}) && isa(comp_model{1}{2}, 'Model')
    comparisonModel = comp_model;
  elseif ~iscell(comp_model)
    comp_model = {comp_model};
  end
  
  % load data
  if isempty(data)
    for d = 1:nData
      dataLoaded = load(dataset{d});
      cmaes_out = getFields(dataLoaded, 'cmaes_out');
      data(d) = cmaes_out.cmaes_out{1}{1};
    end
  end
  
  % load comparison models
  if isempty(comparisonModel)
    for m = 1:nData
      modelLoaded = load(comp_model{m});
      comparisonModel{m} = modelLoaded.models;
    end
    minModels = min(cellfun(@length, comparisonModel));
  end
  
  % load cmaes sampling options
  sOpt = load(fullfile(defDataset, 'bbob_output', 'sampleOpts.mat'));
  sampleOpts = sOpt.sampleOpts;
    
  % testing loop
  for d = 1:nData
    if iscell(dataset)
      fprintf('Testing %s\n', dataset{d})
    else
      fprintf('Testing data %d\n', d)
    end
    actualComp = comparisonModel{d};
    
    % find models with lowest training likelihood
    compLik = cellfun(@(x) x.trainLikelihood, actualComp);
    [~, bestLikId] = sort(compLik);
    gen = bestLikId(1:nGen);
    
    % generation loop
    for g = 1:nGen
      fprintf('Generation: %d\n', gen(g))
      compareModel(g) = actualComp{gen(g)};
      stateVar = compareModel(g).stateVariables;
      trainModel(g) = OrdGpModel(modelOpts, stateVar.xmean');
      % transform coordinates
      if compareModel(g).transformCoordinates
        X = (stateVar.sigma*stateVar.BD*compareModel(g).dataset.X')';
      else
        X = compareModel(g).dataset.X;
      end
      % model training
      trainModel(g) = trainModel(g).train(...
                           X, ...
                           compareModel(g).dataset.y, ...
                           compareModel(g).stateVariables, ...
                           sampleOpts);
      % data for prediction
      genId = data(d).generationStarts(gen(g) + 1);
      Xpred = data(d).arxvalids(:, genId : genId+stateVar.lambda - 1)';
      % model predictions
      yComp{g} = compareModel(g).predict(Xpred);
      yModel{g} = trainModel(g).predict(Xpred);
      
      % get model predictions ordering
      [~, yCompOrd{g}] = sort(yComp{g});
      [~, yModelOrd{g}] = sort(yModel{g});
      
      nTrain = length(yCompOrd{g});
      errDiff(g).mae = sum(abs(yCompOrd{g} - yModelOrd{g})) / nTrain;
      errDiff(g).mzoe = sum(yCompOrd{g} ~= yModelOrd{g}) / nTrain;
      fprintf('GP mae: %f   GP mzoe: %f\n', errDiff(g).mae, errDiff(g).mzoe)
    end
    
    % model likelihood
    modelLik = arrayfun(@(x) trainModel(x).getTrainLikelihood, 1:length(trainModel));
    
    % save models and its predictions
    outputFile = fullfile(outputFolder, ['model_', modelHash(trainModel), '_data',  num2str(d)]);
    if saveIter
      iter = length(dir([outputFile, '*.mat']))+1;
      outputFile = [outputFile, '_', num2str(iter)];
    end
    save(outputFile, 'compareModel', 'trainModel', 'yComp', 'yModel', 'dataset', 'errDiff', 'modelLik')
    % plot results
%     plotLik(compLik(gen), modelLik, gen);

    if d==1 && (sum(modelLik) < Inf)
      overallLik = 0;
    elseif isinf(sum(modelLik))
      return
    end
    % sum of likelihoods
    overallLik = overallLik + sum(modelLik);
  end
  
  fprintf('Sum of model likelihood: %f\n', overallLik)

end

function f = plotLik(compLik, modelLik, gen)
% plot likelihood of the compared models

  ids = 1:length(gen);

  f = figure();
  hold on
  % plot reference model
  plot(ids, compLik, 'g')
  % plot tested model
  plot(ids, modelLik, 'r')
  
  set(gca, 'XTick', ids)
  set(gca, 'XTickLabels', cellstr(num2str(gen'))')
  xlabel('Generation')
  ylabel('Likelihood')
  legend('GpModel', 'OrdGpModel')
  hold off
  
end

function hash = modelHash(m)
% function creating hash for model identification

  % maximum number of characters
  maxChar = 3;

  clm = class(m);
  % model class
  clm = clm(1:maxChar);
  % prediction type
  predT = m.getOptions.prediction(1);
  % binning type
  binT = m.getOptions.binning([1,4]);
  % number of bins
  nob = m(1).getOptions.nBins;
  if ischar(nob)
    nob = strrep(nob, ' ', '');
    nob = nob(1:maxChar);
  else
    nob = num2str(nob);
  end
  
  % create hash
  hash = [clm, predT, binT, nob];
end

function modelOpts = modelOptsParse(mOptVector)
% mOpts = modelOptsParse(mOptVector) parses vector input to model options
% structure

  % TODO: input check
  
  
    % prediction
  switch mOptVector(1)
    case 1
      modelOpts.prediction = 'metric';
    case 2
      modelOpts.prediction = 'avgord';
    otherwise
      return
  end
  % binning
  switch mOptVector(2)
    case 1
      modelOpts.binning = 'none';
    case 2
      modelOpts.binning = 'logcluster';
    case 3
      modelOpts.binning = 'quantile';
    otherwise
      return
  end
  % number of bins
  modelOpts.nBins = ['ceil(', num2str(mOptVector(3)), '*lambda)'];
  % sigma2 hyperparameter
  modelOpts.hyp.lik = mOptVector(4);
  modelOpts.likBounds = mOptVector(5:6);
  % cov hyperparameters
  modelOpts.hyp.cov = mOptVector(7:8);
  modelOpts.covBounds = [mOptVector(9:10); mOptVector(11:12)];
  
end

function sf = getFields(s, varargin)
% sf = getFields(s, fields) extracts fields from structure s

  sf = struct();
  
  if isstruct(s)
    for f = 1:length(varargin)
      if isfield(s, varargin{f})
        sf.(varargin{f}) = s.(varargin{f});
      else
        sf.(varargin{f}) = [];
      end
    end
  end
end