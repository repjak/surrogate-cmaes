function compareModels(modelFolders, func, dims)
% compareModels(modelFolders, func, dims) compares models saved in 
% 'modelFolders' in functions 'func' and dimensions 'dims'.

  % path settings
  defFolder = fullfile('exp', 'experiments', 'model');
  defDatasetPath = fullfile(defFolder, 'defData', 'defSet.mat');
  defModelFolder = fullfile(defFolder, 'defData', 'defModel');
 
  if nargin < 3
    if nargin < 2
      if nargin < 1
        help compareModels
        return
      end
      func = 1;
    end
    dims = 2;
  end
  
  % input check
  if ~iscell(modelFolders)
    modelFolders = {modelFolders};
  end
  nModel = length(modelFolders);
  assert(isnumeric(func), '''func'' has to be integer')
  assert(isnumeric(dims), '''dims'' has to be integer')
  for m = 1:nModel
    assert(isdir(modelFolders{m}), '%s is not folder', modelFolders{m})
  end
  
  nFun = length(func);
  nDim = length(dims);
  
  mse     = cell(nModel, nFun, nDim);
  kendall = cell(nModel, nFun, nDim);
  rde     = cell(nModel, nFun, nDim);
  model   = cell(nModel, nFun, nDim);
  ym      = cell(nModel, nFun, nDim);
  
  % load results
  for m = 1:nModel
    for f = 1:nFun
      for d = 1:nDim
        [~, modelName] = fileparts(modelFolders{m});
        isFile = dir(sprintf('%s*_f%d_%dD.mat', [modelFolders{m}, filesep, modelName], func(f), dims(d)));
        % multiple file case
        if length(isFile) > 1
          warning('Multiple files for function %d dimension %d model %s. Using first one.', func(f), dims(d), modelName)
          isFile = isFile(1);
        end
        % file exists
        if ~isempty(isFile)
          S = load(fullfile(modelFolders{m}, isFile.name));
          model{m, f, d} = S.model;
          ym{m, f, d} = S.ym;
          
          mse_tmp = NaN(size(S.mse));
          kendall_tmp = NaN(size(S.kendall));
          rde_tmp = NaN(size(S.rde));
          % only model not considered as constant will replace NaN with statistics
          isTrained = (arrayfun(@(x) S.model{x}.trainGeneration, 1:length(S.model)) > 0);
          mse_tmp(isTrained) = S.mse(isTrained);
          kendall_tmp(isTrained) = S.kendall(isTrained);
          rde_tmp(isTrained) = S.rde(isTrained);
          
          mse{m, f, d} = mse_tmp;
          kendall{m, f, d} = kendall_tmp;
          rde{m, f, d} = rde_tmp;
        end
      end
    end
  end
  
  % compute statistics
  meanrde = cellfun(@nanmean, rde);
  stdrde  = cellfun(@nanstd, rde);
  meankendall = cellfun(@nanmean, kendall);
  stdkendall  = cellfun(@nanstd, kendall);
  
  fprintf('*** RDE ***\n')
  dispResults(meanrde, func, dims)
  fprintf('*** Kendall ***\n')
  dispResults(meankendall, func, dims)
  
  %TODO: create rankDiff/kendall table

end

function dispResults(stat, func, dims)
% display statistics in table
  for d = 1:length(dims)
    fprintf('  %dD  ', dims(d))
    fprintf('%s\n', strjoin(arrayfun(@(x) sprintf('model_%02d', x), ...
                            1:size(stat,1), 'UniformOutput', false)))
    for f = 1:length(func)
      fprintf(' f%02d ', func(f))
      for m = 1:size(stat, 1)
        if stat(m, f, d) < 0
          fprintf('  %0.3f ', stat(m, f, d))
        else
          fprintf('%s%0.3f ', ones(1, 3-max(0, floor(log10(stat(m, f, d)))))*' ', stat(m, f, d))
        end
      end
      fprintf('\n')
    end
    fprintf('\n')
  end
end