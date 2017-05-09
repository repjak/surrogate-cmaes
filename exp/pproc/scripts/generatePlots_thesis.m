%% ITAT 2016 article plots
% Script for making graphs showing the dependence of minimal function
% values on the number of function evaluations.

%% load data

% checkout file containing all loaded data
tmpFName = fullfile('/tmp', 'thesis_data.mat');
if (exist(tmpFName', 'file'))
  load(tmpFName);
else
  
% default COCO pp colors
coco_colors = struct(...
  'navy', [0, 0, 128], ...
  'magenta', [255, 0, 255], ...
  'orange', [255, 165, 0], ...
  'cornflowerblue', [100, 149, 237], ...
  'red', [255, 0, 0], ...
  'yellowgreen', [154, 205, 50], ...
  'cyan', [0, 255, 255], ...
  'grey', [190, 190, 190] ...
);

% needed function and dimension settings
funcSet.BBfunc = 1:24;
funcSet.dims = [2, 3, 5, 10, 20];

funcSetS.BBfunc = 1:24;
funcSetS.dims = [2, 3, 5, 10, 20];

% folder for results
actualFolder = pwd;
thesisFolder = '/home/jakub/Documents/dp';
plotResultsFolder = fullfile(thesisFolder, 'img/plots');
tableFolder = fullfile(thesisFolder, 'tables');

% path settings
exppath = fullfile('exp', 'experiments');

maes_1_10_path = fullfile(exppath, 'exp_maesEC_15_1_10_cmaes_20D');
%maes_2_10_path = fullfile(exppath, 'exp_maesEC_14_2_10_cmaes_20D');
maes_2_10_path = fullfile(exppath, 'exp_maesEC_16_2_10_cmaes_20D');
maes_3_35_path = fullfile(exppath, 'exp_maesEC_13_5_35_cmaes_20D');
cmaes_path = fullfile(exppath, 'IPOP-CMA-ES-2010');
gpop_path = fullfile(exppath, 'exp_gpop_15_10D_grouped_by_instances');
gpop_10D_path = fullfile(exppath, 'exp_gpop_16_10D');
gpop_20D_path = fullfile(exppath, 'exp_gpop_16_20D');

% scmaes_gp1_path = fullfile(exppath, 'S-CMA-ES_GP1');
% scmaes_gp5_path = fullfile(exppath, 'S-CMA-ES_GP5');
scmaes_10D_path = fullfile(exppath, 'exp_geneEC_10');
scmaes_20D_path = fullfile(exppath, 'exp_geneEC_10_20D');
% scmaes_ada_kl_10D_path = fullfile(exppath, 'exp_geneEC_adapt_14_kl_10D');
% scmaes_ada_kl_20D_path = fullfile(exppath, 'exp_geneEC_adapt_14_kl_20D');
% scmaes_ada_kndl_10D_path = fullfile(exppath, 'exp_geneEC_adapt_15_kendall_10D');
% scmaes_ada_kndl_20D_path = fullfile(exppath, 'exp_geneEC_adapt_15_kendall_20D');
% scmaes_ada_rd_10D_path = fullfile(exppath, 'exp_geneEC_adapt_16_rankdiff_10D');
% scmaes_ada_rd_20D_path = fullfile(exppath, 'exp_geneEC_adapt_16_rankdiff_20D');
% scmaes_ada_rd_f18_20D_path = fullfile(exppath, 'exp_geneEC_adapt_16_rankdiff_f18_20D');
scmaes_ada_kl_10D_path = fullfile(exppath, 'exp_geneEC_adapt_17_kl_10D');
scmaes_ada_kl_20D_path = fullfile(exppath, 'exp_geneEC_adapt_17_kl_20D');
%scmaes_ada_kl_sym_10D_path = fullfile(exppath, 'exp_geneEC_adapt_17_kl_sym_10D');
%scmaes_ada_kl_sym_20D_path = fullfile(exppath, 'exp_geneEC_adapt_17_kl_sym_20D');
scmaes_ada_kndl_10D_path = fullfile(exppath, 'exp_geneEC_adapt_18_kendall_10D');
scmaes_ada_kndl_20D_path = fullfile(exppath, 'exp_geneEC_adapt_18_kendall_20D');
scmaes_ada_rd_10D_path = fullfile(exppath, 'exp_geneEC_adapt_19_rankdiff_10D');
scmaes_ada_rd_20D_path = fullfile(exppath, 'exp_geneEC_adapt_19_rankdiff_20D');

% load data for traditional algorithms
dataFolders = {cmaes_path; ...
               maes_2_10_path; ...
               gpop_10D_path; ...
               gpop_20D_path};

[evals, settings, exp_results] = catEvalSet(dataFolders, funcSet);

dataFolders = {cmaes_path; ...
               scmaes_10D_path; ...
               scmaes_20D_path; ...
               scmaes_ada_kl_10D_path; ...
               scmaes_ada_kl_20D_path; ...
               scmaes_ada_kndl_10D_path; ...
               scmaes_ada_kndl_20D_path; ...
               scmaes_ada_rd_10D_path; ...
               scmaes_ada_rd_20D_path};%cl; ...
               %scmaes_ada_rd_f18_20D_path};

% load data for adaptive scmaes and scmaes
[evals_scmaes, settings_scmaes, exp_results_scmaes] = catEvalSet(dataFolders, funcSetS);

% find ids in settings
clear findSet
findSet.evoControl = 'maes';
findSet.modelOpts.predictionType = 'poi';
findSet.PopSize = 10;
maes_2_10_poi_Id = getStructIndex(settings, findSet);

clear findSet
findSet.evoControl = 'maes';
findSet.modelOpts.predictionType = 'poi';
findSet.PopSize = 35;
maes_5_35_poi_Id = getStructIndex(settings, findSet);

clear findSet
findSet.evoControl = 'maes';
findSet.modelOpts.predictionType = 'fvalues';
findSet.PopSize = 10;
maes_2_10_mmp_Id = getStructIndex(settings, findSet);

clear findSet
findSet.evoControl = 'maes';
findSet.modelOpts.predictionType = 'fvalues';
findSet.PopSize = 35;
maes_5_35_mmp_Id = getStructIndex(settings, findSet);

clear findSet
findSet.gpop_funHistLen = 5;
gpop_Id = getStructIndex(settings, findSet);

clear findSet
findSet.algName = 'IPOP-CMA-ES-2010';
cma_Id = getStructIndex(settings, findSet);

clear findSet
findSet.algName = 'IPOP-CMA-ES-2010';
cma_Id2 = getStructIndex(settings_scmaes, findSet);

clear findSet
%findSet.algName = 'S-CMA-ES_GP1';
findSet.evoControl = 'generation';
findSet.modelType = 'gp';
findSet.evoControlModelGenerations = 1;
scmaes_gp1_Id = getStructIndex(settings_scmaes, findSet);

clear findSet
%findSet.algName = 'S-CMA-ES_GP5';
findSet.evoControl = 'generation';
findSet.modelType = 'gp';
findSet.evoControlModelGenerations = 5;
scmaes_gp5_Id = getStructIndex(settings_scmaes, findSet);

clear findSet
findSet.updaterType = 'KullbackLeibler';
findSet.geneECAdaptive_maxModelGenerations = 5;
scmaes_ada_kl_Id = getStructIndex(settings_scmaes, findSet);

clear findSet
findSet.updaterType = 'KullbackLeibler';
findSet.geneECAdaptive_maxModelGenerations = 10;
scmaes_ada_kl_sym_Id = getStructIndex(settings_scmaes, findSet);

clear findSet
findSet.updaterType = 'rmseKendall';
scmaes_ada_kndl_Id = getStructIndex(settings_scmaes, findSet);

clear findSet
findSet.updaterType = 'RankDiff';
scmaes_ada_rd_Id = getStructIndex(settings_scmaes, findSet);

% extract data
maes_2_10_poi_data = evals(:, :, maes_2_10_poi_Id);
maes_5_35_poi_data = evals(:, :, maes_5_35_poi_Id);
maes_2_10_mmp_data = evals(:, :, maes_2_10_mmp_Id);
maes_5_35_mmp_data = evals(:, :, maes_5_35_mmp_Id);
gpop_data = evals(:, :, gpop_Id);
cmaes_data = evals(:, :, cma_Id);

cmaes_data2 = evals_scmaes(:, :, cma_Id2);
scmaes_gp1_data = evals_scmaes(:, :, scmaes_gp1_Id);
scmaes_gp5_data = evals_scmaes(:, :, scmaes_gp5_Id);
scmaes_ada_kl_data = evals_scmaes(:, :, scmaes_ada_kl_Id);
%scmaes_ada_kl_sym_data = evals_scmaes(:, :, scmaes_ada_kl_sym_Id);
scmaes_ada_kndl_data = evals_scmaes(:, :, scmaes_ada_kndl_Id);
scmaes_ada_rd_data = evals_scmaes(:, :, scmaes_ada_rd_Id);

% extract experiment info such as timing
maes_2_10_poi_info = exp_results(:, :, maes_2_10_poi_Id);
maes_5_35_poi_info = exp_results(:, :, maes_5_35_poi_Id);
maes_2_10_mmp_info = exp_results(:, :, maes_2_10_mmp_Id);
maes_5_35_mmp_info = exp_results(:, :, maes_5_35_mmp_Id);
gpop_info = exp_results(:, :, gpop_Id);

scmaes_gp1_info = exp_results_scmaes(:, :, scmaes_gp1_Id);
scmaes_gp5_info = exp_results_scmaes(:, :, scmaes_gp5_Id);
scmaes_ada_kl_info = exp_results_scmaes(:, :, scmaes_ada_kl_Id);
%scmaes_ada_kl_sym_info = exp_results_scmaes(:, :, scmaes_ada_kl_sym_Id);
scmaes_ada_kndl_info = exp_results_scmaes(:, :, scmaes_ada_kndl_Id);
scmaes_ada_rd_info = exp_results_scmaes(:, :, scmaes_ada_rd_Id);

if (~exist(tmpFName, 'file'))
  save(tmpFName);
end

end

%% Traditional algorithm comparison: CMA-ES, MA-ES-MMP, MA-ES-POI

datanames = {'CMA-ES', 'MA-ES-MMP', 'MA-ES-POI', 'GPOP'};

data = {cmaes_data, maes_2_10_mmp_data, maes_2_10_poi_data, gpop_data};
info_trad = {maes_2_10_mmp_info, maes_2_10_poi_info, gpop_info};

% default COCO pp colors
% colors = [ getAlgColors('cmaes'); [ ....
%   [255 0 255]; ... %[236 0 140]; ...  
%   [255 165 0]; ... %[245 129 55]; ...  
%   [100 149 237]]]; %[65 156 228]] / 255;
colors = [coco_colors.navy; coco_colors.magenta; ...
  coco_colors.orange; coco_colors.cornflowerblue];

maxEvals = 250;

plotFuns = 1:24;
plotDims = [2 3 5 10 20];

for plotDim = plotDims
  clear pdfNames
  pdfNames = {};
  for f = plotFuns
    for d = plotDim
      pdfNames{end+1} = fullfile(plotResultsFolder, sprintf('plot_trad_f%d_%dD', f, d));
    end
  end

  close all
  han = relativeFValuesPlot(data, ...
    'DataNames', datanames, 'DataDims', funcSet.dims, ...
    'DataFuns', funcSet.BBfunc, 'Colors', colors, ...
    'PlotFuns', plotFuns, 'PlotDims', plotDim, ...
    'AggregateDims', false, ...
    'Statistic', @median, 'AggregateFuns', false, ...
    'MaxEval', maxEvals, ...
    'LegendOption', 'first', ...
    'OneFigure', false, ...
    'OmitYLabel', true, ...
    'FunctionNames', true);
  
  print2pdf(han, pdfNames, 1);
end

%% Aggregated traditional algorithms comparison
% Aggregated  scaled function values in dimensions 2, 5, 10.

clear pdfNames
pdfNames = {};
for d = plotDims
  pdfNames{end+1} = fullfile(plotResultsFolder, sprintf('plot_trad_%dD', d));
end

close all
han = relativeFValuesPlot(data, ...
                              'DataNames', datanames, 'DataDims', funcSet.dims, ...
                              'DataFuns', funcSet.BBfunc, 'Colors', colors, ...
                              'PlotFuns', plotFuns, 'PlotDims', plotDims, ...
                              'AggregateDims', false, 'OneFigure', false, ...
                              'Statistic', @median, 'AggregateFuns', true, ...
                              'LineSpecification', {'-', '-', '-', '-', '-', '-', '-'}, ...
                              'LegendOption', 'first', 'MaxEval', maxEvals, ...
                              'OmitYLabel', true, ...
                              'FunctionNames', true);

print2pdf(han, pdfNames, 1)

%% Adaptivity comparison

data_ada = {cmaes_data2, scmaes_gp1_data, scmaes_gp5_data, ...
            scmaes_ada_kl_data, scmaes_ada_kndl_data, scmaes_ada_rd_data};
info_ada = {scmaes_ada_kl_info, scmaes_ada_kndl_info, scmaes_ada_rd_info};

datanames = {'CMA-ES', 'GP-1', 'GP-5', 'ADA-KL', 'ADA-Kendall', 'ADA-RD'};

% colors = [ ....
%   getAlgColors('cmaes'); ...
%   getAlgColors(2); ...
%   getAlgColors('scmaes'); ...
%   getAlgColors(3); ...
%   getAlgColors(4); ...
%   getAlgColors(11); ...
% ];

colors = [coco_colors.navy; coco_colors.magenta; ...
  coco_colors.orange; coco_colors.cornflowerblue; ...
  coco_colors.red; coco_colors.yellowgreen];

maxEvals = 250;

plotFuns = funcSetS.BBfunc;
plotDims = [2 3 5 10 20];

for plotDim = plotDims
  clear pdfNames
  pdfNames = {};
  for f = plotFuns
    for d = plotDim
      pdfNames{end+1} = fullfile(plotResultsFolder, sprintf('plot_ada2_f%d_%dD', f, d));
    end
  end

  close all
  han = relativeFValuesPlot(data_ada, ...
    'DataNames', datanames, 'DataDims', funcSetS.dims, ...
    'DataFuns', funcSetS.BBfunc, 'Colors', colors, ...
    'PlotFuns', plotFuns, 'PlotDims', plotDim, ...
    'AggregateDims', false, ...
    'Statistic', @median, 'AggregateFuns', false, ...
    'MaxEval', maxEvals, ...
    'LegendOption', 'first', ...
    'OneFigure', false, ...
    'OmitYLabel', true, ...
    'FunctionNames', true);
  
  print2pdf(han, pdfNames, 1);
end

%% Aggregated adaptive algorithms comparison
% Aggregated  scaled function values in dimensions 2, 5, 10 and 20.

clear pdfNames
pdfNames = {};
for d = plotDims
  pdfNames{end+1} = fullfile(plotResultsFolder, sprintf('plot_ada2_%dD', d));
end

close all
han = relativeFValuesPlot(data_ada, ...
                              'DataNames', datanames, 'DataDims', funcSetS.dims, ...
                              'DataFuns', funcSetS.BBfunc, 'Colors', colors, ...
                              'PlotFuns', plotFuns, 'PlotDims', plotDims, ...
                              'AggregateDims', false, 'OneFigure', false, ...
                              'Statistic', @median, 'AggregateFuns', true, ...
                              'LineSpecification', {'-', '-', '-', '-', '-', '-', '-'}, ...
                              'LegendOption', 'first', 'MaxEval', maxEvals, ...
                              'OmitYLabel', true, ...
                              'FunctionNames', true);

print2pdf(han, pdfNames, 1)

% %% Multiple comparison of traditional algorithms with a statistical posthoc test.
% 
% close all
% 
% datanames = {'CMA-ES', 'MA-ES-MMP', 'MA-ES-POI', 'GPOP'};
% 
% tableFunc = funcSet.BBfunc;
% tableDims = [2, 3, 5, 10, 20];
% 
% resultDuelTable = fullfile(tableFolder, 'tradDuelTable.tex');
% resultStatsTable = fullfile(tableFolder, 'tradStatsTable.tex');
% 
% [table1, ranks1] = duelTable(data, 'DataNames', datanames, ...
%                             'DataFuns', funcSet.BBfunc, 'DataDims', funcSet.dims, ...
%                             'TableFuns', tableFunc, 'TableDims', tableDims, ...
%                             'Evaluations', [1/3 1], ...
%                             'ResultFile', resultDuelTable);
% 
% [stats1, meanRanks1] = multCompStatsTable(data, 'DataNames', datanames, ...
%                        'DataFuns', funcSet.BBfunc, 'DataDims', funcSet.dims, ...
%                        'TableFuns', tableFunc, 'TableDims', tableDims, ...
%                        'Evaluations', [1/3 1], ...
%                        'ResultFile', resultStatsTable);
% 
% %% Multiple comparison of adaptivity with a statistical posthoc test.
% 
% close all
% 
% datanames = {'CMA-ES', 'GP-1', 'GP-5', 'ADA-KL', 'ADA-Ken', 'ADA-RD'};
% 
% tableFunc = funcSetS.BBfunc;
% tableDims = [2, 3, 5, 10, 20];
% 
% resultDuelTable = fullfile(tableFolder, 'adaDuelTable2.tex');
% resultStatsTable = fullfile(tableFolder, 'adaStatsTable2.tex');
% 
% [table2, ranks2] = duelTable(data_ada, 'DataNames', datanames, ...
%                             'DataFuns', funcSetS.BBfunc, 'DataDims', funcSetS.dims, ...
%                             'TableFuns', tableFunc, 'TableDims', tableDims, ...
%                             'Evaluations', [1/3 1], ...
%                             'ResultFile', resultDuelTable);
% 
% [stats2, meanRanks2] = multCompStatsTable(data_ada, 'DataNames', datanames, ...
%                        'DataFuns', funcSetS.BBfunc, 'DataDims', funcSetS.dims, ...
%                        'TableFuns', tableFunc, 'TableDims', tableDims, ...
%                        'Evaluations', [1/3 1], ...
%                        'ResultFile', resultStatsTable);


% %% CPU timing
% timing_names_trad = {'MAES-MMP', 'MAES-POI', 'GPOP'};
% [~, timing_trad] = cpuTimingTable(info_trad, ...
%                    'DataNames', timing_names_trad, ...
%                    'DataFuns', funcSet.BBfunc, 'DataDims', funcSet.dims);
% 
% for alg = 1:length(info_trad)
%   fprintf(1, '%s ', timing_names_trad{alg});
%   fprintf(1, ' & $%.4f$', mean(timing_trad(:, :, alg)));
%   fprintf(1, '\\\\\n');
% end
% 
% timing_names_ada = {'ADA-KL', 'ADA-Kendall', 'ADA-RD'};
% [~, timing_ada] = cpuTimingTable(info_ada, ...
%                    'DataNames', timing_names_trad, ...
%                    'DataFuns', funcSet.BBfunc, 'DataDims', funcSet.dims);
% 
% for alg = 1:length(info_trad)
%   fprintf(1, '%s ', timing_names_ada{alg});
%   fprintf(1, ' & $%.4f$', mean(timing_ada(:, :, alg)));
%   fprintf(1, '\\\\\n');
% end

% for alg = 1:length(info_trad)
%   cpu_timing = zeros(length(funcSet.BBfunc), length(funcSet.dims));
%   for f = 1:length(funcSet.BBfunc)
%     for d = 1:length(funcSet.dims)
%       info = info_trad{alg};
%       info = info{f, d}{:};
%       nevals = sum(arrayfun(@(x) size(x{:}, 1), info.y_evals));
%       cpu_timing(f, d) = info.time / nevals;
%     end
%   end
%   mean(cpu_timing)
% end
% 
% for alg = 1:length(info_ada)
%   cpu_timing = zeros(length(funcSetS.BBfunc), length(funcSetS.dims));
%   for f = 1:length(funcSetS.BBfunc)
%     for d = 1:length(funcSetS.dims)
%       info = info_ada{alg};
%       nevals = 0;
%       t = 0;
%       for i = 1:length(info{f,d})
%         results = info{f,d}{i};
%         nevals = nevals + sum(arrayfun(@(x) size(x{:}, 1), results.y_evals));
%         t = t + results.time;
%       end
%       cpu_timing(f, d) = t / nevals;
%     end
%   end
%   fprintf(1, ' & $%.4f$', mean(cpu_timing));
%   fprintf(1, '\\\\\n');
% end

                     
%% final clearing
close all
