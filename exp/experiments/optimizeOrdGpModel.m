function [minx, miny] = optimizeOrdGpModel(optimType)
% ordGpModel settings testing script

  defDataFolder = fullfile('exp', 'experiments', 'modelComparison');
  
  lData   = load(fullfile(defDataFolder, 'defData.mat'), 'data');
  lModels = load(fullfile(defDataFolder, 'defModels.mat'), 'comparisonModel');

  dataset = lData.data;
  comp_model = lModels.comparisonModel;
  
%   dataset = {fullfile(defDataset, 'exp_doubleEC_21_short_log_results_2_5D_1.mat'), ...
%              fullfile(defDataset, 'exp_doubleEC_21_short_log_results_3_5D_2.mat'), ...
%              fullfile(defDataset, 'exp_doubleEC_21_short_log_results_6_5D_3.mat'), ...
%              fullfile(defDataset, 'exp_doubleEC_21_short_log_results_8_5D_4.mat'), ...
%              fullfile(defDataset, 'exp_doubleEC_21_short_log_results_13_5D_5.mat'), ...
%              fullfile(defDataset, 'exp_doubleEC_21_short_log_results_17_5D_6.mat'), ...
%              fullfile(defDataset, 'exp_doubleEC_21_short_log_results_21_5D_7.mat')};
% 
%   comp_model = {fullfile(defDataset, 'bbob_output', 'exp_doubleEC_21_short_log_modellog_2_5D_1.mat'), ...
%                 fullfile(defDataset, 'bbob_output', 'exp_doubleEC_21_short_log_modellog_3_5D_2.mat'), ...
%                 fullfile(defDataset, 'bbob_output', 'exp_doubleEC_21_short_log_modellog_6_5D_3.mat'), ...
%                 fullfile(defDataset, 'bbob_output', 'exp_doubleEC_21_short_log_modellog_8_5D_4.mat'), ...
%                 fullfile(defDataset, 'bbob_output', 'exp_doubleEC_21_short_log_modellog_13_5D_5.mat'), ...
%                 fullfile(defDataset, 'bbob_output', 'exp_doubleEC_21_short_log_modellog_17_5D_6.mat'), ...
%                 fullfile(defDataset, 'bbob_output', 'exp_doubleEC_21_short_log_modellog_21_5D_7.mat')};
               
  % default values
  minx = [];
  miny = [];

  % x = [predict, binType, nBinsFactor, lik, liklb, likub, cov1, cov2, cov1lb, cov1ub, cov2lb, cov2ub];
  % inequality constraints
  [Aineq, bineq] = createIneqMatrix({[5, 4, 6], [9, 7, 10], [11, 8, 12]}, 12);
  % lower bounds
  lb = [1, 2, 0, log(1e-10)*ones(1, 9)];
  % upper bounds
  ub = [2, 3, 5, log(1e4)*ones(1, 9)];
  % starting point
  x0 = [2, 2, 1, log(0.01), log(1e-6), log(1e2), log(0.5), log(2), -2, 2, -2, 2];
      
  switch optimType
    case 'fmincon'
      useId = (3:12);
      preset = [2, 1];
      % remove first two columns
      outputFolder =  fullfile('exp', 'experiments', 'modelComparison', ['fmincon_set_', num2str(preset, '%x'), 'x']);

      % x = [nBinsFactor, lik, liklb, likub, cov1, cov2, cov1lb, cov1ub, cov2lb, cov2ub];
      modOptFcn = @(x) compareModels([preset, x], dataset, comp_model, outputFolder, true);

      optimopts = optimoptions( ...
        @fmincon, ...
        'GradObj', 'off', ...
        'Display', 'off', ...
        'MaxIter', 1e2, ...
        'Algorithm', 'interior-point' ...
      );
      optproblem = struct( ...
                  'solver', 'fmincon', ...
                  'objective', modOptFcn, ...
                  'x0', x0(useId), ...
                  'lb', lb(useId), ...
                  'ub', ub(useId), ...
                  'A', Aineq(:, useId), ...
                  'b', bineq, ...
                  'options', optimopts ...
                );

      fprintf(2, 'Existing results of the previous optimization process will be overwritten!\n');
      fprintf('Ensure to store those if necessary.\n')
      [minx, miny, exitflag, optinfo] = fmincon(optproblem);

      % save results
      modOptVect = [preset, minx];
      save(fullfile(outputFolder, 'optim_res'), 'modOptVect', 'minx', 'miny', 'exitflag', 'optinfo')
    % genetic algorithm
    case 'ga'
      preset = [];
      outputFolder =  fullfile('exp', 'experiments', 'modelComparison', ['ga_set_', num2str(preset, '%x'), 'x']);

      % x = [predict, binType, nBinsFactor, lik, liklb, likub, cov1, cov2, cov1lb, cov1ub, cov2lb, cov2ub];
      modOptFcn = @(x) compareModels([preset, x], dataset, comp_model, outputFolder, true);
      
      gaopt = gaoptimset(@ga);
      gaopt = gaoptimset(gaopt, 'PopulationSize', 1, 'InitialPopulation', x0);
      
      [minx, miny, exitflag, output] = ...
        ga(modOptFcn, 12, Aineq, bineq, [], [], lb, ub, [], [1, 2], gaopt);

      % save results
      modOptVect = [preset, minx];
      save(fullfile(outputFolder, 'optim_res'), 'modOptVect', 'minx', 'miny', 'exitflag', 'output')
    case 'cmaes'
    otherwise
      return
  end
end

function [A, b] = createIneqMatrix(ineqVect, nVar)
% creates inequality matrix from simplified syntax:
%  [2,5,4] -> x2 < x5 < x4 
%
% Input:
%   ineqVect - cell-array of double vectors
%   nVar     - number of variables in inequalities

  if nargin < 2
    nVar = 0;
  end
  
  nVar = max(max(max(cell2mat(ineqVect))), nVar);
  nEq = cellfun(@length,ineqVect)-1;
  nSet = length(nEq);
  sumEq = sum(nEq);
  
  A = zeros(sumEq, nVar);
  b = zeros(sumEq, 1);
  
  e = 0;
  for s = 1:nSet
    for ine = 1:nEq(s)
      e = e + 1;
      A(e, ineqVect{s}(ine))   =  1;
      A(e, ineqVect{s}(ine+1)) = -1;
    end
  end
end