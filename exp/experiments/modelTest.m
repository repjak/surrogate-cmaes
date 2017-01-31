function [mse, kendall, rde, model, ym] = modelTest(modelType, modelOpts, ds)
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

  mse = zeros(nDatasetsPerInstance, 1);
  kendall = zeros(nDatasetsPerInstance, 1);
  rde = zeros(nDatasetsPerInstance, 1);
  model = cell(nDatasetsPerInstance, 1);
  ym = cell(nDatasetsPerInstance, 1);

  for i = 1:nDatasetsPerInstance
    m = ModelFactory.createModel(modelType, modelOpts, ds.means{i});
    m = m.train(ds.trainSetX{i}, ds.trainSetY{i}, ds.cmaesStates{i}, ds.sampleOpts);

    ym{i} = m.predict(ds.testSetX{i});

    mse(i) = sum((ym{i} - ds.testSetY{i}).^2) / length(ym{i});
    kendall(i) = corr(ym{i}, ds.testSetY{i}, 'type', 'kendall');
    rde(i) = errRankMu(ym{i}, ds.testSetY{i}, floor(length(ym{i})/2));
    
    model{i} = m;

    fprintf('Model (gen. # %3d) MSE = %e, Kendall = %.2f, rankDiffErr = %.2f\n', ...
      ds.generations(i), mse(i), kendall(i), rde(i));
  end
end