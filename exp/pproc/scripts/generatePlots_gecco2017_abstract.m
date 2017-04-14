%% GECCO 2017 abstract plots and tables
% Script for making graphs showing the dependence of minimal function
% values on the number of function values and tables showing differences
% between ORDGP and GP DTS-CMA-ES.
% 
% Created for GECCO 2017 abstract.

%% initial settings

tmpFName = fullfile('/tmp', 'gecco2017abstract_data.mat');
if (exist(tmpFName', 'file'))
  load(tmpFName);
else

% settings
func = (1:24);
dims = [2, 5, 10];
maxEvals = 100;
% needed function and dimension settings for GRAPHS !!!
funcSet.BBfunc = func;
funcSet.dims = [2, 5];

% folder for results
actualFolder = pwd;
articleFolder = fullfile(actualFolder(1:end - 1 - length('surrogate-cmaes')), 'latex_scmaes', 'gecco2017abstract');
plotResultsFolder = fullfile(articleFolder, 'images');
tableFolder = fullfile(articleFolder, 'tex');
if ~isdir(plotResultsFolder)
  mkdir(plotResultsFolder)
end
if ~isdir(tableFolder)
  mkdir(tableFolder)
end

% path settings
exppath = fullfile('exp', 'experiments');
defModelFolder = fullfile(exppath, 'model');

ord_2D_path = fullfile(exppath, 'exp_doubleEC_ord_05_2D');
ord_2D_none_path = fullfile(exppath, 'exp_doubleEC_ord_05_2D_none');
ord_5D_path = fullfile(exppath, 'exp_doubleEC_ord_04_5D');
ord_5D_none_path = fullfile(exppath, 'exp_doubleEC_ord_05_5D_none');

cmaes_path = fullfile(exppath, 'CMA-ES');
saacmes_path = fullfile(exppath, 'BIPOP-saACM-k');
dts_path = fullfile(exppath, 'DTS-CMA-ES_05_2pop');
lmm_path = fullfile(exppath, 'lmm-CMA-ES');

% load data from online testing
dataFolders = {ord_2D_path; ...
               ord_2D_none_path; ...
               ord_5D_path; ...
               ord_5D_none_path; ...
               cmaes_path; ...
               saacmes_path; ...
               dts_path; ...
               lmm_path};
             
[evals, settings] = catEvalSet(dataFolders, funcSet);

% find ids in settings
clear findSet
findSet.modelOpts.binning = 'none';
se_noneId = getStructIndex(settings, findSet);

findSet.modelOpts.binning = 'logcluster';
findSet.modelOpts.nBins = 'mu';
se_clust_muId = getStructIndex(settings, findSet);
findSet.modelOpts.nBins = 'lambda';
se_clust_lamId = getStructIndex(settings, findSet);

findSet.modelOpts.binning = 'unipoints';
findSet.modelOpts.nBins = 'mu';
se_unip_muId = getStructIndex(settings, findSet);
findSet.modelOpts.nBins = 'lambda';
se_unip_lamId = getStructIndex(settings, findSet);

clear findSet
findSet.algName = 'CMA-ES';
cma_Id = getStructIndex(settings, findSet);
findSet.algName = 'BIPOP-saACM-k';
saacm_Id = getStructIndex(settings, findSet);
findSet.algName = 'DTS-CMA-ES_05_2pop';
dts_Id = getStructIndex(settings, findSet);
findSet.algName = 'lmm-CMA-ES';
lmm_Id = getStructIndex(settings, findSet);

% extract data
se_av_non_data = evals(:, :, se_noneId);
se_av_clu_mu_data = evals(:, :, se_clust_muId);
se_av_clu_lam_data = evals(:, :, se_clust_lamId);
se_av_uni_mu_data = evals(:, :, se_unip_muId);
se_av_uni_lam_data = evals(:, :, se_unip_lamId);
cmaes_data = evals(:, :, cma_Id);
saacmes_data = evals(:, :, saacm_Id);
dtscmaes_data = evals(:, :, dts_Id);
lmmcmaes_data = evals(:, :, lmm_Id);

% color settings
seAvNonCol    = getAlgColors(1);
seAvCluMuCol  = getAlgColors(2);
seAvCluLamCol = getAlgColors(3);
seAvUniMuCol  = getAlgColors(4);
seAvUniLamCol = getAlgColors(5);

cmaesCol   = getAlgColors('cmaes');
saacmesCol = getAlgColors('saacmes');
dtsCol     = getAlgColors('dtscmaes');
lmmCol     = getAlgColors('lmmcmaes');

if (~exist(tmpFName, 'file'))
  save(tmpFName);
end

end
             
%% Algorithm comparison: ordGP, DTS-CMA-ES, lmm-CMA-ES, saACMES, CMA-ES
% Aggregation of function values across dimensions 2, 5.

data = {se_av_uni_lam_data, ...
        lmmcmaes_data, ...
        saacmes_data, ...
        dtscmaes_data, ...
        cmaes_data};

datanames = {'Ord-DTS', 'lmm-CMA-ES', 'BIPOP-{}^{s*}ACMES-k', 'DTS-CMA-ES', 'CMA-ES'};

colors = [seAvUniLamCol; lmmCol; saacmesCol; dtsCol; cmaesCol]/255;

plotDims = [2, 5];

clear pdfNames
pdfNames = {fullfile(plotResultsFolder, 'alg2_5D')};

close all
han = relativeFValuesPlot(data, ...
                              'DataNames', datanames, 'DataDims', funcSet.dims, ...
                              'DataFuns', funcSet.BBfunc, 'Colors', colors, ...
                              'PlotFuns', funcSet.BBfunc, 'PlotDims', plotDims, ...
                              'AggregateDims', false, 'OneFigure', true, ...
                              'Statistic', @median, 'AggregateFuns', true, ...
                              'LineSpecification', {'-.', '-.', '-', '-', '-', '-'}, ...
                              'LegendOption', 'split', 'MaxEval', 100);

                              
                            
print2pdf(han, pdfNames, 1)

%% Algorithm comparison on chosen functions: ordGP, DTS-CMA-ES, lmm-CMA-ES,
% saACMES, CMA-ES
% Aggregation of function values across dimensions 2, 5.

data = {se_av_uni_lam_data, ...
        lmmcmaes_data, ...
        saacmes_data, ...
        dtscmaes_data, ...
        cmaes_data};

datanames = {'Ord-DTS', 'lmm-CMA-ES', 'BIPOP-{}^{s*}ACMES-k', 'DTS-CMA-ES', 'CMA-ES'};

colors = [seAvUniLamCol; lmmCol; saacmesCol; dtsCol; cmaesCol]/255;

plotDims = [2, 5];
plotFuns = [6];

clear pdfNames
pdfNames = fullfile(plotResultsFolder, 'alg_f6_2_5D');

close all
han = relativeFValuesPlot(data, ...
                              'DataNames', datanames, 'DataDims', funcSet.dims, ...
                              'DataFuns', funcSet.BBfunc, 'Colors', colors, ...
                              'PlotFuns', plotFuns, 'PlotDims', plotDims, ...
                              'AggregateDims', false, 'OneFigure', true, ...
                              'Statistic', @median, 'AggregateFuns', false, ...
                              'LineSpecification', {'-.', '-.', '-', '-', '-', '-'}, ...
                              'LegendOption', 'hide', 'MaxEval', 100, ...
                              'FunctionNames', true);

                              
                            
print2pdf(han, pdfNames, 1)

%% zip resulting plot files

zip([plotResultsFolder, '.zip'], fullfile(plotResultsFolder, '*.pdf'))