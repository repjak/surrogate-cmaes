classdef GenerationsUpdaterFactory
  % Creates a GenerationsUpdater object based on the parameter
  % updaterParams.updaterType

  methods (Static)
    function obj = createUpdater(ec, surrogateOpts)
      if ~isfield(surrogateOpts, 'updaterType') || isempty(surrogateOpts.updaterType)
       obj = GenerationsUpdaterFactory.createDefaultUpdater(surrogateOpts);
       return;
      end

      switch lower(surrogateOpts.updaterType)
        case 'rmsekendall'
          obj = GenerationsUpdaterRMSEKendall(ec, surrogateOpts);
        case 'rankdiff'
          obj = GenerationsUpdaterRankDiff(ec, surrogateOpts);
        case 'kullbackleibler'
          obj = GenerationsUpdaterKL(ec, surrogateOpts);
        otherwise
          % including surrogateOpts.updaterType == 'constant'
          obj = GenerationsUpdaterFactory.createDefaultUpdater(ec, surrogateOpts);
      end
    end

    function obj = createDefaultUpdater(ec, surrogateOpts)
      % a fallback into a constant updater, i.e. the generation-based EC
      % without adaptation
      if (~isfield(surrogateOpts, 'geneECAdaptive_modelGenerations') || ...
          isempty(surrogateOpts.geneECAdaptive_modelGenerations) || ...
          ~isfield(surrogateOpts, 'geneECAdaptive_origGenerations') || ...
          isempty(surrogateOpts.geneECAdaptive_origGenerations))
        if ~(isfield(surrogateOpts, 'evoControlOrigGenerations') && ...
            isfield(surrogateOpts, 'evoControlModelGenerations'))
          error('There''s not enough parameters for GenerationsUpdaterConstant');
        else
          p = struct('origGenerations', surrogateOpts.evoControlOrigGenerations, ...
                'modelGenerations', surrogateOpts.evoControlModelGenerations);
        end
      else
        p = surrogateOpts;
      end

      obj = GenerationsUpdaterConstant(ec, p);
    end
  end
end
