%% GECCO 2017 BBOB workshop plots and tables
% Script for making graphs showing the dependence of minimal function
% values on the number of function values and tables showing differences
% between ORDGP and GP DTS-CMA-ES.
% 
% Created for GECCO 2017 BBOB workshop article.

%% initial settings

tmpFName = fullfile('/tmp', 'gecco2017bbob_workshop_data.mat');
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
articleFolder = fullfile(actualFolder(1:end - 1 - length('surrogate-cmaes')), 'latex_scmaes', 'gecco2017bbob');
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

% model names
% {covMaterniso, 5}
mat_avg_none      = 'ordgpModel_814382';
% mat_met_none      = 'ordgpModel_814470';
mat_avg_clus_mu   = 'ordgpModel_982019';
mat_avg_clus_lam  = 'ordgpModel_1026675';
mat_avg_clus_2lam = 'ordgpModel_1034955';
mat_avg_unip_mu   = 'ordgpModel_971627';
mat_avg_unip_lam  = 'ordgpModel_1016028';
mat_avg_unip_2lam = 'ordgpModel_1024280';
% mat_met_clus_mu   = 'ordgpModel_982107';
% mat_met_clus_lam  = 'ordgpModel_1026763';
% mat_met_clus_2lam = 'ordgpModel_1035043';
% mat_met_unip_mu   = 'ordgpModel_971715';
% mat_met_unip_lam  = 'ordgpModel_1016116';
% mat_met_unip_2lam = 'ordgpModel_1024368';
% squaredexponential
se_avg_none      = 'ordgpModel_585187';
se_avg_clus_mu   = 'ordgpModel_727377';
se_avg_clus_lam  = 'ordgpModel_765672';
se_avg_clus_2lam = 'ordgpModel_770672';
se_avg_unip_mu   = 'ordgpModel_718580';
se_avg_unip_lam  = 'ordgpModel_756620';
se_avg_unip_2lam = 'ordgpModel_761592';

% gp model folder
gpModelFolder = fullfile(defModelFolder, 'defData', ['defModel_', num2str(maxEvals), 'FE']);
% modelFNames = {mat_avg_none; mat_met_none; ...
%                 mat_avg_clus_mu; mat_avg_clus_lam; mat_avg_clus_2lam; ...
%                 mat_avg_unip_mu; mat_avg_unip_lam; mat_avg_unip_2lam; ...
%                 mat_met_clus_mu; mat_met_clus_lam; mat_met_clus_2lam; ...
%                 mat_met_unip_mu; mat_met_unip_lam; mat_met_unip_2lam};
modelFileNames = {mat_avg_none; ... 
                mat_avg_clus_mu; mat_avg_clus_lam; mat_avg_clus_2lam; ...
                mat_avg_unip_mu; mat_avg_unip_lam; mat_avg_unip_2lam; ...
               se_avg_none; ...
                se_avg_clus_mu; se_avg_clus_lam; se_avg_clus_2lam; ...
                se_avg_unip_mu; se_avg_unip_lam; se_avg_unip_2lam};
              
% ordgp model folders
modelFolders = cellfun(@(x) fullfile(defModelFolder, [x, '_', num2str(maxEvals), 'FE']), ...
                       modelFileNames, 'UniformOutput', false);
modelFolders = [gpModelFolder; modelFolders];

              
% compute model statistics
stats = modelStatistics(modelFolders, func, dims, false);

% load data from online testing
dataFolders = {ord_2D_path; ...
               ord_2D_none_path; ...
               ord_5D_path; ...
               ord_5D_none_path; ...
               cmaes_path; ...
               dts_path}; % ...
%               saacmes_path; ...
%               lmm_path};
             
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
% findSet.algName = 'BIPOP-saACM-k';
% saacm_Id = getStructIndex(settings, findSet);
findSet.algName = 'DTS-CMA-ES_05_2pop';
dts_Id = getStructIndex(settings, findSet);
% findSet.algName = 'lmm-CMA-ES';
% lmm_Id = getStructIndex(settings, findSet);

% extract data
se_av_non_data = evals(:, :, se_noneId);
se_av_clu_mu_data = evals(:, :, se_clust_muId);
se_av_clu_lam_data = evals(:, :, se_clust_lamId);
se_av_uni_mu_data = evals(:, :, se_unip_muId);
se_av_uni_lam_data = evals(:, :, se_unip_lamId);
cmaes_data = evals(:, :, cma_Id);
% saacmes_data = evals(:, :, saacm_Id);
dtscmaes_data = evals(:, :, dts_Id);
% lmmcmaes_data = evals(:, :, lmm_Id);

% color settings
seAvNonCol    = getAlgColors(1);
seAvCluMuCol  = getAlgColors(2);
seAvCluLamCol = getAlgColors(3);
seAvUniMuCol  = getAlgColors(4);
seAvUniLamCol = getAlgColors(5);

cmaesCol   = getAlgColors('cmaes');
% saacmesCol = getAlgColors('saacmes');
dtsCol     = getAlgColors('dtscmaes');
% lmmCol     = getAlgColors('lmmcmaes');

if (~exist(tmpFName, 'file'))
  save(tmpFName);
end

end

%% create statistic tables
% model names
gpdts  = '\dtsgp';
covMat = '\covMat';
covSE  = '\covSE';
clus   = '\setHAC';
unip   = '\setQuant';
mu     = '\mu';
lam    = '\lambda';
lam2   = '2\lambda';
modelNames = {gpdts; ...
              covMat; ... 
              strjoin({covMat, clus, mu}, ', '); ...
              strjoin({covMat, clus, lam}, ', '); ...
              strjoin({covMat, clus, lam2}, ', '); ...
              strjoin({covMat, unip, mu}, ', '); ...
              strjoin({covMat, unip, lam}, ', '); ...
              strjoin({covMat, unip, lam2}, ', '); ...
              covSE; ... 
              strjoin({covSE, clus, mu}, ', '); ...
              strjoin({covSE, clus, lam}, ', '); ...
              strjoin({covSE, clus, lam2}, ', '); ...
              strjoin({covSE, unip, mu}, ', '); ...
              strjoin({covSE, unip, lam}, ', '); ...
              strjoin({covSE, unip, lam2}, ', ')};
% add $ to model names
modelNames = cellfun(@(x) ['$ ', x, ' $'], modelNames, 'UniformOutput', false);

modelStatTable(stats, ...
               'DataDims', [2, 5, 10],...
               'TableDims', 5, ...
               'Format', 'tex', ...
               'ResultFolder', tableFolder, ...
               'ModelNames', modelNames, ...
               'ShowCaption', false, ...
               'ShowDimStat', true)