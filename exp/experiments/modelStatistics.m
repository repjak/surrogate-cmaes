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

  % multicomparison results
  alpha = 0.05;
  stats.multcomprde = multcomp(rde, nModel, alpha);
  stats.multcompkendall = multcomp(kendall, nModel, alpha);
  stats.multcomprankmzoe = multcomp(rankmzoe, nModel, alpha);

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

function [b, validreps] = cell2blockmat(stat, reps)
% Prepare a "block" matrix for Friedman's test. Blocked variables are
% functions and dimensions, i.e. the columns of the resulting matrix
% correspond to model settings and rows correspond to repetitions over
% all functions and dimensions.
%
% Input:
%   stat - a performance result
%   reps - number of repetitions for each combination of blocking variables
%
% Output:
%   b         - a matrix with repetitions in rows blocked by functions
%               and dimensions
%   validreps - minimum number of non-NaN repetitions over all combinations
%               of blocking variables
  if nargin < 2
    reps = length(stat{1, 1});
  end

  b = zeros(size(stat, 2) * size(stat, 3) * reps, size(stat, 1));
  validreps = reps;

  for m = 1:size(stat, 1)
    for d = 1:size(stat, 3)
      for f = 1:size(stat, 2)
        s = stat{m, f, d};
        validreps = min(validreps, sum(~isnan(s)));
        for i = 1:reps
          j = i;
          while isnan(s(j))
            j = j+1;
          end
          b((f-1)*reps + (d-1)*f*reps + j, m) = s(i);
        end
      end
    end
  end
end

function a = multcomp(res, n, alpha)
% Compute multicompare statistics for a peformance metrics.
%
% Input:
%   res   - a performance metrics results
%   n     - a number of models
%   alpha - a significance level
%
% Output:
%   a     - a logical 1-n array with 1s on the indices of those models that
%           significantly differ from all other (n-1) models
  [~, reps] = cell2blockmat(res);
  % the second pass is to remove NaNs while keeping consistent number of
  % repetitions for blocking variables
  b = cell2blockmat(res, reps);
  [~, ~, friedstats] = friedman(b, reps, 'off');
  c = multcompare(friedstats, 'CType','bonferroni', 'Display', 'off');
  a = finddiff(c, n, alpha);
  a = a == (n - 1);
end

function a = finddiff(c, n, alpha)
% For each model, find the number of other models which a significant
% difference.
%
% Input:
%   c     - an output from multcompare
%   n     - a number of models
%   alpha - a significance level
%
% Output:
%   a     - an 1xn array with a(i) giving the number of models that differ
%           significantly from ith model
  a = zeros(1, n);
  for i = 1:size(c, 1)
    if c(i, 6) <= alpha
      i1 = int16(c(i, 1));
      i2 = int16(c(i, 2));
      a(i1) = a(i1) + 1;
      a(i2) = a(i2) + 1;
    end
  end
end
