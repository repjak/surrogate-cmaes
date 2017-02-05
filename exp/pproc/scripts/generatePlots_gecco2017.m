%% GECCO 2017 article plots and tables
% Script for making graphs showing the dependence of minimal function
% values on the number of function values and tables showing differences
% between ORDGP and GP DTS-CMA-ES.
% 
% Created for GECCO 2017 article.

%% initial settings

tmpFName = fullfile('/tmp', 'gecco2017data.mat');
if (exist(tmpFName', 'file'))
  load(tmpFName);
else

% settings
func = (1:24);
dims = [2, 5, 10];
maxEvals = 100;

% folder for results
actualFolder = pwd;
articleFolder = fullfile(actualFolder(1:end - 1 - length('surrogate-cmaes')), 'latex_scmaes', 'gecco2017paper');
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
ord_path = fullfile(exppath, 'exp_doubleEC_ord_02');

% model names
% {covMaterniso, 5}
mat_avg_none      = 'ordgpModel_814382';
mat_met_none      = 'ordgpModel_814470';
mat_avg_clus_mu   = 'ordgpModel_814470';
mat_avg_clus_lam  = 'ordgpModel_1026675';
mat_avg_clus_2lam = 'ordgpModel_1034955';
mat_avg_unip_mu   = 'ordgpModel_971627';
mat_avg_unip_lam  = 'ordgpModel_1016028';
mat_avg_unip_2lam = 'ordgpModel_1024280';
mat_met_clus_mu   = 'ordgpModel_982107';
mat_met_clus_lam  = 'ordgpModel_1026763';
mat_met_clus_2lam = 'ordgpModel_1035043';
mat_met_unip_mu   = 'ordgpModel_971715';
mat_met_unip_lam  = 'ordgpModel_1016116';
mat_met_unip_2lam = 'ordgpModel_1024368';
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

if (~exist(tmpFName, 'file'))
  save(tmpFName);
end

end

%% create statistic tables
% model names
covMat = '\covMat';
covSE  = '\covSE';
clus = '\mathrm{C}';
unip = '\mathrm{Q}';
mu = '\mu';
lam = '\lambda';
lam2 = '2\lambda';
modelNames = {'\mathrm{GP}'; ...
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
               'Format', 'tex', ...
               'ResultFolder', tableFolder, ...
               'ModelNames', modelNames)