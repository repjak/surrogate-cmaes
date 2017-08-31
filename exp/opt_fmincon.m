function [x, y_evals, stopflag, varargout] = opt_fmincon(FUN, DIM, ftarget, maxfunevals, id, varargin)
% minimizes FUN in DIM dimensions by multistarts of fminsearch.
% ftarget and maxfunevals are additional external termination conditions,
% where at most 2 * maxfunevals function evaluations are conducted.
% fminsearch was modified to take as input variable usual_delta to
% generate the first simplex.
% set options, make sure we always terminate
% with restarts up to 2*maxfunevals are allowed
%
% varargin/1 argument can be 'xstart' -- starting point for the optimization
% varargin/2 can be 'cmOptions' -- options for the BOBYQA algorithm


  fDelta = 1e-8;
  stopflag = [];
  varargout = {};

  cmOptions = struct( ...
    'Algorithm', 'interior-point', ...
    ... % 'ObjectiveLimit', ftarget,
    'Display', 'off', ...
    'MaxFunctionEvaluations', min(1e8*DIM, maxfunevals), ...
    'LBounds', -5, ...
    'UBounds',  5 ...
    );

  y_evals = [];

  if (nargin >= 6), exppath = [varargin{1} filesep];
    else              exppath = '';  end
  if (nargin >= 7), xstart = varargin{2};
    else              xstart = 8 * rand(dim, 1) - 4;  end
  % xstart = 8 * rand(DIM, 1) - 4;

  load([exppath 'scmaes_params.mat'], 'bbParamDef', 'sgParamDef', 'cmParamDef');
  [~, ~, cmParams] = getParamsFromIndex(id, bbParamDef, sgParamDef, cmParamDef);

  for fname = fieldnames(cmParams)'
    cmOptions.(fname{1}) = cmParams.(fname{1});
  end

  lb = cmOptions.LBounds * ones(1, DIM);
  ub = cmOptions.UBounds * ones(1, DIM);

  [x, fmin, stopflag, output] = fmincon(FUN, xstart, [], [], [], [], lb, ub, [], cmOptions);
  
  y_evals = [fmin, output.funcCount NaN NaN NaN];

  if (nargout > 3)
    varargout{1} = [];
  else
    varargout = cell(0);
  end
end % function
