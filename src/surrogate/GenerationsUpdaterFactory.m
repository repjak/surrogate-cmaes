classdef GenerationsUpdaterFactory
  %GENERATIOSUPDATERFACTORY Summary of this class goes here
  %   Detailed explanation goes here

  methods (Static)
    function obj = createUpdater(surrogateOpts)
      switch lower(surrogateOpts.updaterParams.updaterType)
        case 'rmse'
          if ~isfield(surrogateOpts, 'updaterParams')
            error('Missing parameters for GenerationsUpdaterRMSE');
          end
          obj = GenerationsUpdaterRMSE(surrogateOpts.updaterParams);
        otherwise
          % including surrogateOpts.updaterType == 'constant'
          %
          % for backward compatibility
          if ~(isfield(surrogateOpts, 'updaterParams')) || isempty(surrogateOpts.updaterParams)
            if ~(isfield(surrogateOpts, 'evoControlOrigGenerations') && ...
                isfield(surrogateOpts, 'evoControlModelGenerations'))
              error('There''s enough parameters for GenerationsUpdaterConstant');
            else
              p = struct('origGenerations', surrogateOpts.evoControlOrigGenerations, ...
                    'modelGenerations', surrogateOpts.evoControlModelGenerations);
            end
          else
            p = surrogateOpts.updaterParams;
          end

          obj = GenerationsUpdaterConstant(p);
      end
    end
  end
end