function ds = modelTestSets(exp_id, fun, dim)
% ds = modelTestSets(exp_id, fun, dim) creates/loads dataset from 
% experiment
%
% Input:
%   exp_id - experiment id
%   fun    - functions to load
%   dim    - dimensions to load
%
% Output:
%   ds - loaded data | #fun x #dim struct 

  if nargin < 3
    if nargin < 2
      if nargin < 1
        exp_id = 'exp_doubleEC_21_log';
      end
      fun = 1;
    end
    dim = 2;
  end
  
  nDatasetsPerInstance = 10;

  idStart = 1;

  % path settings
  experimentPath = fullfile('exp', 'experiments', exp_id);
  outputFolder = fullfile('exp', 'experiments', 'model');
  dataFolder = fullfile(outputFolder, 'defData');
  modelFolder = fullfile(dataFolder, 'defModel');
  datasetName = fullfile(dataFolder, 'defSet');
  defModelName = fullfile(modelFolder, 'gpModel');
  
  if ~exist(outputFolder, 'dir')
    mkdir(outputFolder)
  end
  if ~exist(dataFolder, 'dir')
    mkdir(dataFolder)
  end
  if ~exist(modelFolder, 'dir')
    mkdir(modelFolder)
  end

  nFun = length(fun);
  nDim = length(dim);
  
  % dimension loop
  for d = 1:nDim
    % function loop
    for f = 1:nFun
      %TODO: fix id according to dimension
      id = idStart + f - 1;
      fprintf('#### f%d in %dD ####\n', fun(f), dim(d));

      ds(f,d) = datasetFromInstance(exp_id, nDatasetsPerInstance, fun(f), dim(d), id);
      
      modelFile = sprintf('%s_f%d_%dD.mat', defModelName, fun(f), dim(d));
      
      if ~exist(modelFile, 'file')

        scmaesOutFile = sprintf('%s/%s_results_%d_%dD_%d.mat', experimentPath, exp_id, fun(f), dim(d), id);
        load(scmaesOutFile, 'cmaes_out', 'exp_settings', 'surrogateParams');
        
        [mse, kendall, rde, model, ym] = modelTest(surrogateParams.modelType, surrogateParams.modelOpts, ds);
        
        % save model results
        save(modelFile, 'mse', 'kendall', 'rde', 'model', 'ym')
      end
      
    % function loop end  
    end
    
  % dimension loop end
  end
  
  % save default dataset
  save(datasetName, 'ds', 'fun', 'dim')
  
end