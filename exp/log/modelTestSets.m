function ds = modelTestSets(exp_id, fun, dim, maxEval)
% ds = modelTestSets(exp_id, fun, dim, maxEval) creates/loads dataset from 
% experiment.
%
% Input:
%   exp_id  - experiment id
%   fun     - functions to load
%   dim     - dimensions to load
%   maxEval - maximal number of evaluations times dimension to load
%
% Output:
%   ds - loaded data | #fun x #dim cell-array
%
% See Also:
%   modelTest

  if nargin < 4
    if nargin < 3
      if nargin < 2
        if nargin < 1
          exp_id = 'exp_doubleEC_21_log';
        end
        fun = 1;
      end
      dim = 2;
    end
    maxEval = 250;
  end

  % default settings
  nDatasetsPerInstance = 10;

  % path settings
  experimentPath = fullfile('exp', 'experiments', exp_id);
  outputFolder = fullfile('exp', 'experiments', 'model');
  dataFolder = fullfile(outputFolder, 'defData');
  datasetName = ['defSet_', num2str(maxEval),'FE'];
  modelFolder = fullfile(dataFolder, ['defModel_', num2str(maxEval),'FE']);
  datasetFile = fullfile(dataFolder, [datasetName '.mat']);
  defModelName = fullfile(modelFolder, 'defModel');

  % create folders
  if ~exist(outputFolder, 'dir')
    mkdir(outputFolder)
  end
  if ~exist(dataFolder, 'dir')
    mkdir(dataFolder)
  end
  if ~exist(modelFolder, 'dir')
    mkdir(modelFolder)
  end

  % check experiment parameters
  if exist(fullfile(experimentPath, 'scmaes_params.mat'), 'file')
    exp_par = load(fullfile(experimentPath, 'scmaes_params.mat'));
    exp_dim = cell2mat(exp_par.bbParamDef(strcmp({exp_par.bbParamDef.name}, 'dimensions')).values);
    exp_fun = cell2mat(exp_par.bbParamDef(strcmp({exp_par.bbParamDef.name}, 'functions')).values);
  else
    warning('No scmaes_params.mat found in %s directory. Using default fun and dim settings.', experimentPath)
    exp_dim = 2;
    exp_fun = 1;
  end

  % check fun and dim values
  inExpFun = ismember(fun, exp_fun);
  inExpDim = ismember(dim, exp_dim);
  if ~all(inExpFun)
    warning('Functions %s are not in experiment %s. Removing from loading.', num2str(fun(~inExpFun), '%d '), exp_id);
    fun = fun(inExpFun);
  end
  if ~all(inExpDim)
    warning('Dimensions %s are not in experiment %s. Removing from loading.', num2str(dim(~inExpDim), '%d '), exp_id);
    dim = dim(inExpDim);
  end

  nFun = length(fun);
  nDim = length(dim);

  f_ds = struct('ds', {}, 'fun', [], 'dim', [], 'maxEval', maxEval);
  if (exist(datasetFile, 'file'))
    fprintf('The dataset file "%s" already existed.\n   Copying to "%s.bak".\n', datasetFile, datasetFile);
    copyfile(datasetFile, [datasetFile '.bak']);
    f_ds = load(datasetFile);
    if (isfield(f_ds, 'ds') && ~isempty(f_ds.ds))
      ds = f_ds.ds;
    end
  else
    ds = cell(nFun, nDim);
  end

  % dimension loop
  for d = 1:nDim
    d_exp = find(dim(d) == exp_dim);
    % function loop
    for f = 1:nFun
      f_exp = find(fun(f) == exp_fun);
      id = (d_exp-1)*length(exp_fun) + f_exp;
      fprintf('#### f%d in %dD ####\n', fun(f), dim(d));

      if (~isempty(ds{f, d}) && isstruct(ds{f, d}))
        fprintf('...already loaded\n');
        continue
      end

      modelFile = sprintf('%s_f%d_%dD.mat', defModelName, fun(f), dim(d));
      if (exist(modelFile, 'file'))
        mfile = load(modelFile); % should be loaded: 'stats', 'model', 'ym', 'ds_actual'
      end
      if (isfield(mfile, 'ds_actual'))
        ds_actual = mfile.ds_actual;
        ds{f, d} = ds_actual;
      else
        % load dataset from saved modellog/cmaes_out of the corresponding instance
        if isempty(ds{f, d})
          ds_actual = datasetFromInstance(exp_id, nDatasetsPerInstance, fun(f), dim(d), id, maxEval);
        end
      end
      % succesfully loaded
      if isstruct(ds_actual)
        ds{f, d} = ds_actual;
        % compute default models if they do not exist
        if ~exist(modelFile, 'file')
          scmaesOutFile = sprintf('%s/%s_results_%d_%dD_%d.mat', experimentPath, exp_id, fun(f), dim(d), id);
          load(scmaesOutFile, 'cmaes_out', 'exp_settings', 'surrogateParams');
          [stats, model, ym] = modelTest(surrogateParams.modelType, surrogateParams.modelOpts, ds{f, d});
          % save model results
          save(modelFile, 'stats', 'model', 'ym', 'ds_actual');
        end
      end
    % function loop end  
    end

  % dimension loop end
  end

  % save default dataset
  fun = union(fun, f_ds.fun);
  dim = union(dim, f_ds.dim);
  maxEval = min(maxEval, f_ds.maxEval);
  save(datasetFile, 'ds', 'fun', 'dim', 'maxEval')

end
