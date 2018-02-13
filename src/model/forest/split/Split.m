classdef Split
% Split creates a split function used in decision trees
  
  properties (Constant)
    splitCandidate = struct( ...
      'splitter', @(X) zeros(size(X, 1), 1), ...
      'gain', -inf, ...
      'leftID', 1, ...
      'rightID', 1 ...
      );
  end
  
  properties %(Access = protected)
    split_transformationOptions % transform options
    split_transformation        % transformation
    split_X                     % input data
    split_y                     % output data
    split_allEqual              % whether all y values are equal
    split_soft                  % use soft split
    split_lambda                % lambda steepness in soft logit function
  end
  
  methods
    function obj = Split(options)
      obj.split_transformationOptions = defopts(options, 'split_transformationOptions', struct);
      obj.split_soft = defopts(options, 'split_soft', false);
      obj.split_lambda = defopts(options, 'split_lambda', 1);
    end
    
    function obj = reset(obj, X, y)
    % sets new transformed input
      [obj.split_X, obj.split_y, obj.split_transformation] = ...
        transform(X, y, obj.split_transformationOptions);
      n = size(X, 1);
      varY = var(y);
      obj.split_allEqual = any(isnan(varY)) || all(varY < eps(max(varY)) * n);
    end
    
    function best = get(obj, splitGain)
    % returns the split with max splitGain
      best = obj.splitCandidate;
      if obj.split_allEqual
        return
      end
      candidate = obj.splitCandidate;
      candidate.gain = splitGain(candidate.splitter);
%       candidate.modelID = 0;
      best = candidate;
    end
    
    function best = getDiscrAnal(obj, splitGain, c, best, discrimTypes)
    % returns split using discriminant analysis
      for i = 1:numel(discrimTypes)
        discrimType = discrimTypes{i};
        try
          model = fitcdiscr(obj.split_X, c, 'DiscrimType', discrimType);
          modelTrained = true;
        catch
          modelTrained = false;
        end
        
        cl = unique(c);
        % singular covariance matrix
        if ~modelTrained && ...
            sum(cl(1) == c) > 1 && sum(cl(2) == c) > 1 && ...
            any(strcmp(discrimType, {'linear', 'quadratic'}))
            
          pseudoDiscrimType = strcat('pseudo', discrimType);
          try
            model = fitcdiscr(obj.split_X, c, 'DiscrimType', pseudoDiscrimType);
            modelTrained = true;
          catch
          end
        end
        
        % model trained => create candidate
        if modelTrained
          candidate.splitter = obj.createModelSplitter(model);
          [candidate.gain, candidate.leftID, candidate.rightID] = splitGain.get(candidate.splitter);
          if candidate.gain > best.gain
            best = candidate;
          end
        end
      end
    end
  end
  
  methods (Access = protected)
    function f = createSplitter(obj, splitter)
      trans = obj.split_transformation;
      if obj.split_soft
        % soft splitter using logit
        % normalize first
        r = splitter(obj.split_X);
        idx = r <= 0;
        rLeft = r(idx);
        rRight = r(~idx);
        mm = [-sum(rLeft) / numel(rLeft), sum(rRight) / numel(rRight)];
        splitter = @(X) Split.normalize(splitter(X), mm);

        lambda = obj.split_lambda;
        f = @(X) 1 ./ (1 + exp(-lambda * splitter(transformApply(X, trans))));
      else
        % hard splitter
        f = @(X) splitter(transformApply(X, trans)) <= 0;
      end
    end
    
    function f = createModelSplitter(obj, model)
      trans = obj.split_transformation;
      if obj.split_soft
        % soft splitter
        f = @(X) Split.modelProbability(model, transformApply(X, trans));
      else
        % hard splitter
        f = @(X) model.predict(transformApply(X, trans)) == 1;
      end
    end
  end
  
  methods (Access = private, Static)
    function p = modelProbability(model, X)
      [~, p] = model.predict(X);
    end
    
    function r = normalize(r, mm)
      idx = r <= 0;
      if mm(1) > 0
        r(idx) = r(idx) / mm(1);
      end
      if mm(2) > 0
        r(~idx) = r(~idx) / mm(2);
      end
    end
  end
end