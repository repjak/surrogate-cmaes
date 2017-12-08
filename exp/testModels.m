function modelFolder = testModels(modelOptions, opts, funcToTest, dimsToTest, instToTest)
% testModels -- tests fitting power of models on DTS dataset.
%
% modelFolder = testModels(...) tests models and returns resulting names 
% of folders.
%
% Input:
%   modelOptions - model options to test | struct (cell-array of struct)
%   opts         - other options and settings
%       .modelType      - type of model (from ECFactory) | string
%       .exp_id         - exp_id of the _model testing_ experiment | string
%       .dataset        - filename (w/o extension) | string
%       .maxEvals       - maximal number of FE per dimensions to consider | integer
%       .exppath_short  - directory with experiments subdirectories | string
%       .statistics     - cell-array of statistics' names | cell array of strings
%       .rewrite_results - whether to rewrite already saved results | bool
%       .saveModels     - whether the models should be saved, too, default FALSE | bool
%       .alwaysRetrain  - let always be the models be retrained (efectively
%                         deletes models in the 'data')
%   funcToTest         - functions to test | array of integers
%   dimsToTest         - dimensions to test | array of integers
%   instToTest         - dimensions to test | array of integers
%
% Output:
%   modelFolder  - list of folders containing results | cell-array of 
%                  string
%
% See Also:
%   modelTestSet, testOneModel

  % Default input parameters settings
  if nargin < 5
    instToTest = 1;
    if nargin < 4
      dimsToTest = 2;
      if nargin < 3
        funcToTest = 1;
      end
    end
  end

  % Ensure input parameters as cell arrays
  modelType = opts.modelType;
  if ~iscell(modelType),        modelType = {modelType}; end
  if ~iscell(modelOptions),     modelOptions = {modelOptions}; end
  nModel = length(modelOptions);
  if length(modelType) == 1 && nModel > 1
    modelType = modelType(ones(nModel, 1));
  end
  assert(isnumeric(funcToTest), '''funcToTest'' has to be integer')
  assert(isnumeric(dimsToTest), '''dimsToTest'' has to be integer')
  assert(isnumeric(instToTest), '''instToTest'' has to be integer')

  % Default options in 'opts'
  opts.maxEvals = defopts(opts, 'maxEvals', 250);
  opts.exppath_short = defopts(opts, 'exppath_short', fullfile('exp', 'experiments'));
  opts.statistics = defopts(opts, 'statistics', { 'mse' });
  opts.saveModels = defopts(opts, 'saveModels', false);

  % Load the dataset and its parameters
  assert(isfield(opts, 'dataset'), 'opts does not include a dataset');
  if ischar(opts.dataset)
    assert(exist(opts.dataset, 'file')==2, 'Dataset %s does not exist', opts.dataset)
    loadData = load(opts.dataset);
  else
    assert(iscell(opts.dataset), 'Dataset has to be string or cell')
    loadData = struct();
  end
  data          = defopts(loadData, 'ds', opts.dataset);
  dataFunc      = defopts(loadData, 'fun', funcToTest);
  dataDims      = defopts(loadData, 'dim', dimsToTest);
  dataInst      = defopts(loadData, 'inst', instToTest);
  dataMaxEvals  = defopts(loadData, 'maxEval', 250);
  dataNSnapshots = defopts(loadData, 'nSnapshotsPerRun', 10);
  maxEvals      = min(dataMaxEvals, opts.maxEvals);

  % filter the required functions/dimensions/instances to those that are really
  % in the dataset
  dimsToTest = restricToDataset(dimsToTest, dataDims, 'Dimensions');
  funcToTest = restricToDataset(funcToTest, dataFunc, 'Functions');
  instToTest = restricToDataset(instToTest, dataInst, 'Instances');

  % Assign (hopefully) unique names and output directories to models
  for m = 1:nModel
    modelHashName{m} = [modelType{m}, 'model_', modelHash(modelOptions{m})];
    % TODO: ensure, that this name is really unique
    modelFolder{m}   = fullfile(opts.exppath_short, opts.exp_id, [modelHashName{m}, '_', num2str(maxEvals), 'FE']);
  end

  % dimension loop
  for dim = dimsToTest
    d_data = find(dim == dataDims, 1);

    % function loop
    for fun = funcToTest
      f_data = find(fun == dataFunc, 1);

      % Check that we really have this data
      if (isempty(data{f_data, d_data}))
        warning('Data of func %d in dim %d is missing in the dataset!', dataFunc(f_data), dataDims(d_data));
        continue;
      end

      % model loop
      for m = 1:nModel
        fprintf('*******************  Fun: %d  Dim: %d  Model: %d  *******************\n', fun, dim, m)

        % directory and filename for the output
        [~, ~] = mkdir(modelFolder{m});
        modelFile = fullfile(modelFolder{m}, sprintf('%s_f%d_%dD.mat', modelHashName{m}, fun, dim));

        % do not rewrite existing files unless wanted
        finishedInstances = false(size(instToTest));
        if (exist(modelFile, 'file') && ~opts.rewrite_results)
          % there are some already calculated results, try to load them and continue
          try
            if (opts.saveModels)
              oldResults = load(modelFile, 'stats', 'models', 'y_models', 'instances', 'modelOptions', 'fun', 'dim');
            else
              oldResults = load(modelFile, 'stats', 'y_models', 'instances', 'modelOptions', 'fun', 'dim');
            end
            finishedInstances = ismember(instToTest, oldResults.instances);
          catch err
            warning('File %s is corrupted. Will be rewritten. Error: %s', modelFile, err.message);
          end
          if (all(finishedInstances))
            fprintf('All instances in %s already calculated. Skipping model testing.\n', modelFile);
            continue;
          end
        end
        if (any(finishedInstances))
          startInstanceIdx = find(~finishedInstances, 1);
          fprintf('Some instances already calculated in %s.\n', modelFile);
          fprintf('Starting from inst. # %d.\n', instToTest(startInstanceIdx));
          stats = oldResults.stats;
          if (opts.saveModels)
            models = oldResults.models;
          else
            models = cell(length(instToTest), dataNSnapshots);
          end
          y_models = oldResults.y_models;
        else
          % (re-)calculate the results
          startInstanceIdx = 1;
          if (exist(modelFile, 'file'))
            warning('Stop testing if you do not want to rewrite file %s', modelFile)
          end

          % prepare output variables
          stats = struct();
          for st = 1:length(opts.statistics)
            stats.(opts.statistics{st}) = NaN(length(instToTest), dataNSnapshots);
          end
          models   = cell(length(instToTest), dataNSnapshots);
          y_models = cell(length(instToTest), dataNSnapshots);
        end

        % print what you are going to test
        printStructure(modelOptions{m});

        % instances loop
        for ii = startInstanceIdx:length(instToTest)
          inst = instToTest(ii);
          i_data = find(inst == dataInst, 1);
          fprintf('-- instance %2d --\n', inst);

          BBOB_DIR = [opts.scratch '/bbob_output'];
          fgeneric('initialize', fun, inst, BBOB_DIR);
          % Info about tested BBOB function
          bbob_handlesF = benchmarks('handles');
          noisyHandles = benchmarksnoisy('handles');
          bbob_handlesF(100+(1:length(noisyHandles))) = noisyHandles;
          opts.bbob_func = bbob_handlesF{fun};

          % delete the saved models from the dataset to be always retrained
          if (isfield(opts, 'alwaysRetrain') && opts.alwaysRetrain)
            if (isfield(data{f_data, d_data, i_data}, 'models'))
              data{f_data, d_data, i_data}.models = cell(size(data{f_data, d_data, i_data}.models));
            end
            if (isfield(data{f_data, d_data, i_data}, 'models2'))
              data{f_data, d_data, i_data}.models2 = cell(size(data{f_data, d_data, i_data}.models2));
            end
          end

          % train & test the model on the 'dataNSnapshots' datasets
          % from the current instance
          if (~isfield(opts, 'trySecondModel') || ~opts.trySecondModel)
            [new_stats, thisModels, y_models(ii, :)] = ...
                testOneModel(modelType{m}, modelOptions{m}, ...
                data{f_data, d_data, i_data}, dataNSnapshots, opts);
            modelsVarString = { 'y_models' };
            if (opts.saveModels)
              modelsVarString(end+1) = 'models';
            end
          else
            [new_stats, thisModels, y_models(ii, :), thisModels2, y_models2(ii, :)] = ...
                testOneModel(modelType{m}, modelOptions{m}, ...
                data{f_data, d_data, i_data}, dataNSnapshots, opts);
            modelsVarString = { 'y_models', 'y_models2' };
            if (opts.saveModels)
              modelsVarString(end+(1:2)) = { 'models', 'models2' };
            end
          end

          % finalize BBOB function settings
          fgeneric('finalize');

          % save results into output variables
          for st = 1:length(opts.statistics)
            fname = opts.statistics{st};
            stats.(fname)(ii, :) = new_stats.(fname);
          end

          % save results of the so-far calculated instances
          instances = instToTest(1:ii);
          if (opts.saveModels)
            models(ii, :) = thisModels;
            if (isfield(opts, 'trySecondModel') && opts.trySecondModel)
              models2(ii, :) = thisModels2;
            end
          end
          save(modelFile, 'stats', modelsVarString{:}, 'instances', 'modelOptions', 'fun', 'dim')
        end  % instance loop
      end  % model loop
    end  % function loop
  end  % dimension loop

end


function [sField, sVal] = getFieldsVals(s)
% sf = getFields(s, fields) extracts fields and its values from structure s

  sField = fieldnames(s);
  nFields = length(sField);
  sVal = cell(nFields, 1);
  for i = 1:nFields
    sVal{i} = s.(sField{i});
  end
end


function [intToTest] = restricToDataset(intToTest, dataInt, identifier)
  id_ok = ismember(intToTest, dataInt);
  if ~all(id_ok)
    warning('%s %s are not in tested dataset', identifier, ...
      strjoin(arrayfun(@num2str, intToTest(~id_ok), 'UniformOutput', false), ','))
    intToTest = intToTest(id_ok);
  end
end
