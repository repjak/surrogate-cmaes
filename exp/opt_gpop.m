function [x, ilaunch, y_evals, stopflag, varargout] = opt_gpop(FUN, dim, ftarget, maxfunevals, id, varargin)
% minimizes FUN in dim dimensions by multistarts of gpop.
% ftarget and maxfunevals are additional external termination conditions,
% where at most 2 * maxfunevals function evaluations are conducted.
%
% fminsearch was modified to take as input variable usual_delta to
% generate the first simplex.
% set options, make sure we always terminate
% with restarts up to 2*maxfunevals are allowed
%
% last varargin argument can be 'exppath' -- path to the experiment's data

% Be aware: 'id' is an additional parameter!

varargout = cell(nargout);

xstart = 8 * rand(dim, 1) - 4; % random start solution

fDelta = 1e-8;

% GPOP defaults
gpopOptions = struct( ...
  'maxFunEvals', min(1e8*dim, maxfunevals), ...
  'stopFitness', ftarget ...
);

% CMA-ES defaults
cmOptions = struct( ...
  'MaxFunEvals', 250*dim, ...
  'LBounds', -5 * ones(size(xstart)), ...
  'UBounds',  5 * ones(size(xstart)), ...
  'LogTime',  0, ...
  'SaveVariables', 'off', ...
  'LogModulo', 0, ...
  'DispModulo', '1000' ...
);

y_evals = [];

if (nargin >= 6)
	exppath = [varargin{1} filesep];
else
	exppath = '';
end

load([exppath 'scmaes_params.mat'], 'bbParamDef', 'sgParamDef', 'cmParamDef', 'exp_id', 'exppath_short', 'logDir');
[bbParams, sgParams, cmParams] = getParamsFromIndex(id, bbParamDef, sgParamDef, cmParamDef);

for fname = fieldnames(cmParams)'
	cmOptions.(fname{1}) = cmParams.(fname{1});
end

% copy params 'gpop_*' from sgParams and trim the prefix from the beginning
prefix = 'gpop_';
l_prefix = length(prefix);
for fname = fieldnames(sgParams)'
  if (length(fname{1}) > l_prefix && ...
      isequal(strfind(fname{1}, prefix), 1))
    name = fname{1};
    gpopOptions.(name(l_prefix+1:end)) = sgParams.(name);
  end
end

% Info about tested function is for debugging purposes
bbob_handlesF = benchmarks('handles');
sgParams.modelOpts.bbob_func = bbob_handlesF{bbParams.functions(1)};
sgParams.expFileID = [num2str(bbParams.functions(1)) '_' num2str(dim) 'D_' num2str(id)];

ilaunch = 1; % no restarts

[x, fmin, counteval, stopflag, y_eval] = gpop(FUN, xstart, gpopOptions, cmOptions, sgParams.modelOpts);

n_y_evals = size(y_eval,1);
y_eval(:,1) = y_eval(:,1) - (ftarget - fDelta) * ones(n_y_evals,1);
y_evals = [y_evals; y_eval];

end % function
