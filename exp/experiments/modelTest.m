function [stats, model, ym] = modelTest(modelType, modelOpts, ds)
% modelTest(modelType, modelOpts, ds) compute statistics on chosen model
%
% Input:
%   modelType - type of tested model | string
%   modelOpts - options of tested model | struct
%   ds        - dataset for testing | struct
%
% Output:
%   mse     - mean square error
%   kendall - Kendall tau rank
%   rde     - ranking difference error
%   m       - trained models | cell-array
%   ym      - predicted values | cell-array
%
% See Also:
%   datasetFromInstance

  nDatasetsPerInstance = length(ds.trainSetX);

  mse = NaN(nDatasetsPerInstance, 1);
  mzoe = NaN(nDatasetsPerInstance, 1);
  kendall = NaN(nDatasetsPerInstance, 1);
  rde = NaN(nDatasetsPerInstance, 1);
  model = cell(nDatasetsPerInstance, 1);
  ym = cell(nDatasetsPerInstance, 1);

  for i = 1:nDatasetsPerInstance
    m = ModelFactory.createModel(modelType, modelOpts, ds.means{i});
    m = m.train(ds.trainSetX{i}, ds.trainSetY{i}, ds.cmaesStates{i}, ds.sampleOpts);

    if m.isTrained()
      y = ds.testSetY{i};
      ym{i} = m.predict(ds.testSetX{i});
      % calculate statistics
      mse(i)     = predictionStats(y, ym{i}, 'mse');
      mzoe(i)    = predictionStats(y, ym{i}, 'mzoe');
      kendall(i) = predictionStats(y, ym{i}, 'kendall');
      rde(i)     = predictionStats(y, ym{i}, 'rde');
    end
    
    model{i} = m;

    fprintf('Model (gen. # %3d) MSE = %e, Kendall = %.2f, rankDiffErr = %.2f\n', ...
      ds.generations(i), mse(i), kendall(i), rde(i));
  end
  
  % create stats structure
  stats.mse = mse;
  stats.mzoe = mzoe;
  stats.kendall = kendall;
  stats.rde = rde;
end