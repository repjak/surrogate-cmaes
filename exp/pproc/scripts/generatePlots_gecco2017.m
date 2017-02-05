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
dims = [2, 5, 10, 20];
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
avg_none      = 'ordgpModel_814382';
met_none      = 'ordgpModel_814470';
avg_clus_mu   = 'ordgpModel_814470';
avg_clus_lam  = 'ordgpModel_1026675';
avg_clus_2lam = 'ordgpModel_1034955';
avg_unip_mu   = 'ordgpModel_971627';
avg_unip_lam  = 'ordgpModel_1016028';
avg_unip_2lam = 'ordgpModel_1024280';
met_clus_mu   = 'ordgpModel_982107';
met_clus_lam  = 'ordgpModel_1026763';
met_clus_2lam = 'ordgpModel_1035043';
met_unip_mu   = 'ordgpModel_971715';
met_unip_lam  = 'ordgpModel_1016116';
met_unip_2lam = 'ordgpModel_1024368';
% model folders
gpModelFolder = fullfile(defModelFolder, 'defData', ['defModel_', num2str(maxEvals), 'FE']);
modelNames = {avg_none; met_none; ...
                avg_clus_mu; avg_clus_lam; avg_clus_2lam; ...
                avg_unip_mu; avg_unip_lam; avg_unip_2lam; ...
                met_clus_mu; met_clus_lam; met_clus_2lam; ...
                met_unip_mu; met_unip_lam; met_unip_2lam};
modelFolders = cellfun(@(x) fullfile(defModelFolder, [x, '_', num2str(maxEvals), 'FE']), ...
                       modelNames, 'UniformOutput', false);
modelFolders = [gpModelFolder; modelFolders];
              
% compute model statistics
stats = modelStatistics(modelFolders, func, dims, false);

if (~exist(tmpFName, 'file'))
  save(tmpFName);
end

end

%% compare results
compareModels(modelFolders, func, dims)