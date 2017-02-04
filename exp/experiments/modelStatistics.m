function stats = modelStatistics(modelFolders, func, dims, dispRes)
% stats = modelStatistics(modelFolders, func, dims, dispRes) computes model
% statistics Kendall rank and ranking difference error.
%
% Input:
%   modelFolders - folders containing model results
%   func         - functions for statistics computation
%   dims         - dimensions for statistics computation
%   dispRes      - display results | boolean
 
  if nargin < 4
    if nargin < 3
      if nargin < 2
        if nargin < 1
          help modelStatistics
          return
        end
        func = 1;
      end
      dims = 2;
    end
    dispRes = true;
  end
  
  % input check
  if ~iscell(modelFolders)
    modelFolders = {modelFolders};
  end
  nModel = length(modelFolders);
  assert(isnumeric(func), '''func'' has to be integer')
  assert(isnumeric(dims), '''dims'' has to be integer')
  for m = 1:nModel
    assert(isdir(modelFolders{m}), '%s is not a folder', modelFolders{m})
  end
  
  nFun = length(func);
  nDim = length(dims);
  
  mse      = cell(nModel, nFun, nDim);
  mzoe     = cell(nModel, nFun, nDim);
  kendall  = cell(nModel, nFun, nDim);
  rankmse  = cell(nModel, nFun, nDim);
  rankmzoe = cell(nModel, nFun, nDim);
  rde      = cell(nModel, nFun, nDim);
  model    = cell(nModel, nFun, nDim);
  ym       = cell(nModel, nFun, nDim);
  
  % load results
  for m = 1:nModel
    for f = 1:nFun
      for d = 1:nDim
        [~, modelName] = fileparts(modelFolders{m});
        isFile = dir(sprintf('%s*_f%d_%dD.mat', [modelFolders{m}, filesep, modelName], func(f), dims(d)));
        % modelName can end with _#FE
        if isempty(isFile)
          modelNameEnd = regexp(modelName, '._\d*FE');
          if ~isempty(modelNameEnd)
            isFile = dir(sprintf('%s*_f%d_%dD.mat', ...
              [modelFolders{m}, filesep, modelName(1:modelNameEnd)], ...
              func(f), dims(d)));
          end
        end
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
          
          mse_tmp      = NaN(size(S.stats.mse));
          mzoe_tmp     = NaN(size(S.stats.mzoe));
          kendall_tmp  = NaN(size(S.stats.kendall));
          rankmse_tmp  = NaN(size(S.stats.rankmse));
          rankmzoe_tmp = NaN(size(S.stats.rankmzoe));
          rde_tmp      = NaN(size(S.stats.rde));
          % only model not considered as constant will replace NaN with statistics
          isTrained = (arrayfun(@(x) S.model{x}.trainGeneration, 1:length(S.model)) > 0);
          mse_tmp(isTrained)      = S.stats.mse(isTrained);
          mzoe_tmp(isTrained)     = S.stats.mzoe(isTrained);
          kendall_tmp(isTrained)  = S.stats.kendall(isTrained);
          rankmse_tmp(isTrained)  = S.stats.rankmse(isTrained);
          rankmzoe_tmp(isTrained) = S.stats.rankmzoe(isTrained);
          rde_tmp(isTrained)      = S.stats.rde(isTrained);
          
          mse{m, f, d}      = mse_tmp;
          mzoe{m, f, d}     = mzoe_tmp;
          kendall{m, f, d}  = kendall_tmp;
          rankmse{m, f, d}  = rankmse_tmp;
          rankmzoe{m, f, d} = rankmzoe_tmp;
          rde{m, f, d}      = rde_tmp;
        end
      end
    end
  end
  
  % compute statistics
  stats.meanrde      = cellfun(@nanmean, rde);
  stats.stdrde       = cellfun(@nanstd,  rde);
  stats.meankendall  = cellfun(@nanmean, kendall);
  stats.stdkendall   = cellfun(@nanstd,  kendall);
  stats.meanrankmzoe = cellfun(@nanmean, rankmzoe);
  stats.stdrankmzoe  = cellfun(@nanstd,  rankmzoe);
  
  % display results
  if dispRes
    fprintf('\n')
    fprintf('*** RDE ***\n')
    dispResults(stats.meanrde, func, dims)
    fprintf('*** Kendall ***\n')
    dispResults(stats.meankendall, func, dims)
    fprintf('*** MZOE ***\n')
    dispResults(stats.meanrankmzoe, func, dims)
    
    % if no output is required do not return anything
    if nargout == 0
      clear stats
    end
  end

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
        if isnan(stat(m, f, d))
          fprintf('    ---  ')
        elseif stat(m, f, d) < 0
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