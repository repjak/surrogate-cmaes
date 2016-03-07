exp_id = 'exp_gpop_10D';
exp_description = 'GPOP (Buche) with GP and CMA-ES, 24, functions, 15 instances';

login = 'repjak';
logMatlabOutput = true;
logDir = '/storage/plzen1/home/repjak/public';

% BBOB parameters
% which dimensions to optimize, subset of [2 3 5 10 20 40];
bbParamDef(1).name   = 'dimensions';
bbParamDef(1).values = {2, 5, 10};      % {2, 5 10};
% function ID's to optimize (2 Sphere, 3 Rastrigin, 8 Rosenbrock)
bbParamDef(2).name   = 'functions';
bbParamDef(2).values = num2cell(1:24);  % {1, 2, 3, 5, 6, 8, 10, 11, 12, 13, 14, 20, 21};
% function being optimized -- BBOB wrap-around with header
% xbest = function( fun, dim, ftarget, maxfunevals )
bbParamDef(3).name   = 'opt_function';
bbParamDef(3).values = {@opt_gpop};
bbParamDef(4).name   = 'instances';
bbParamDef(4).values = {[1:5 31:40]};   % default is [1:5, 31:40]
bbParamDef(5).name   = 'maxfunevals';   % MAXFUNEVALS - 10*dim is a short test-experiment
bbParamDef(5).values = {'250 * dim'};   % increment maxfunevals successively

% GP model options
sgParamDef(1).name   = 'modelOpts';
% hyperparameters
% ell = 0.5 * ones(dim, 1); % characteristic length scale
% sf2 = 0.1;                % signal variance
% sf = 0.01;                % constant shift
% sn = 1e-4;                % signal noise
sgParamDef(1).values = { struct( ...
  'trainAlgorithm', 'fmincon', ...
  'covFcn', '{@covSum, {@covSEard, @covConst, @covNoise}}', ...
  'hyp', struct('lik', log(0.001), ...
    'cov', 'log([0.5 * ones(dim, 1); 0.1; 0.01; 1e-4])'), ...
  'normalizeY', true) };

% GPOP parameters, historically part of sgParam
% TODO: allow defining new parameter groups in experiment definition
% number of training points selected by distance to xbest
sgParamDef(2).name   = 'gpop_nc';
sgParamDef(2).values = {'5 * dim'};
% number of training points selected by time of evaluation
sgParamDef(3).name   = 'gpop_nr';
sgParamDef(3).values = {'5 * dim'};
% parameters of merit function (added density measure)
sgParamDef(4).name   = 'gpop_meritParams';
sgParamDef(4).values = {[0 1 2 4]};
% add perturbance if solution doesn't improve more than tolXPert
sgParamDef(5).name   = 'gpop_tolXPrtb';
sgParamDef(5).values = {1e-3};

% CMA-ES parameters
cmParamDef(1).name   = 'PopSize';
cmParamDef(1).values = {'(4 + floor(3*log(N)))' '10'};
cmParamDef(2).name   = 'ParentNumber';
cmParamDef(2).values = {2};
cmParamDef(3).name   = 'Restarts';
cmParamDef(3).values = {4};

% path to current file -- do not change this
pathstr = fileparts(mfilename('fullpath'));
exppath  = [pathstr filesep exp_id];
exppath_short  = pathstr;
[s,mess,messid] = mkdir(exppath);
[s,mess,messid] = mkdir([exppath filesep 'cmaes_results']);
addpath(exppath);

% Save the directory of the experiment data for debugging purposes
sgParamDef(end+1).name = 'experimentPath';
sgParamDef(end).values = { exppath };

% TODO: allow saving to 'gpop_params.mat' 
save([exppath filesep 'scmaes_params.mat'], 'bbParamDef', 'sgParamDef', 'cmParamDef', 'exp_id', 'exppath_short', 'logDir');

% run the rest of the scripts generation
generateShellScriptsMetacentrum
