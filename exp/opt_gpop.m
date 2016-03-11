function [x, ilaunch, y_evals, stopflag, varargout] = opt_gpop(FUN, dim, ftarget, maxfunevals, id, varargin)
% minimizes FUN in DIM dimensions by multistarts of fminsearch.
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

xstart = 8 * rand(dim, 1) - 4; % random start solution

% fDelta = 1e-8;

% GPOP defaults
gpopOptions = struct( ...
  'stopMaxFunEvals', min(1e8*dim, maxfunevals), ...
  'stopMaxIter', 40 * maxfunevals, ...
  'stopMaxIterPrtb', 10);

% CMA-ES defaults
% TODO: termination criteria for individual CMA-ES runs
%  'MaxFunEvals', min(1e8*DIM, maxfunevals), ...
%  'StopFitness', ftarget, ...
cmOptions = struct( ...
  'MaxFunEvals', min(10*dim, maxfunevals), ...
  'LBounds', -5, ...
  'UBounds',  5, ...
  'LogTime',  0, ...
  'SaveVariables', 'off', ...
  'LogModulo', 0, ...
  'DispModulo', '10');

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

for ilaunch = 1:1e4
  % Info about tested function is for debugging purposes
  bbob_handlesF = benchmarks('handles');
  sgParams.modelOpts.bbob_func = bbob_handlesF{bbParams.functions(1)};
  sgParams.expFileID = [num2str(bbParams.functions(1)) '_' num2str(dim) 'D_' num2str(id)];

  [x, fmin, counteval, stopflag] = gpop(FUN, xstart, gpopOptions, cmOptions, sgParams.modelOpts);

  % TODO: y_eval
  varargout = cell(0);

  if (feval(FUN, 'fbest') < ftarget || ...
      feval(FUN, 'evaluations') >= maxfunevals)
    break;
  end

  xstart = x;
end % for

  function stop = callback(x, optimValues, state)
    stop = false;
    if optimValues.fval < ftarget
      stop = true;
    end
  end % function callback

end % function