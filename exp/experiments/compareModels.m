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
        isFile = dir(sprintf('%s*_f%d_%dD.mat', modelFolders{m}, func(f), dims(d)));
        % multiple file case
        if length(isFile) > 1
          warning('Multiple files for function %d dimension %d. Using first one.', func(f), dims(d))
          isFile = isFile(1);
        end
        % file exists
        if ~isempty(isFile)
          S = load(fullfile(modelFolders{m}, isFile.name));
          mse{m, f, d} = S.mse;
          kendall{m, f, d} = S.kendall;
          rde{m, f, d} = S.rde;
          model{m, f, d} = S.model;
          ym{m, f, d} = S.ym;
        end
      end
    end
  end
  
  %TODO: create rankDiff/kendall table

end