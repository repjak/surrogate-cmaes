classdef PairObliqueSplit < Split
% PairObliqueSplit for each pair of points constructs a normal vector,
% projects all points onto this vector and uses the projected values as
% tresholds for decision to which side the point belongs

  properties %(Access = protected)
    split_nQuantize % quantization of tresholds
  end

  methods
    function obj = PairObliqueSplit(options)
      obj = obj@Split(options);
      obj.split_nQuantize = defopts(options, 'split_nQuantize', 0);
    end
    
    function best = get(obj, splitGain)
    % returns the split with max splitGain
      best = obj.splitCandidate;
      if obj.split_allEqual
        return
      end
      [n, ~] = size(obj.split_X);
      for i = 1:n-1
        for j = i+1:n
          % normal vector of the hyperplane
          v = obj.split_X(i, :) - obj.split_X(j, :);
          % project X onto v
          values = (obj.split_X * v')';
          if obj.split_nQuantize > 0 && numel(values) > obj.split_nQuantize
            mm = minmax(values);
            tresholds = linspace(mm(1), mm(2), obj.split_nQuantize);
          else
            tresholds = unique(values);
          end
          for treshold = tresholds
            candidate = obj.splitCandidate;
            candidate.splitter = obj.createSplitter(@(X) ...
              X * v' - treshold);
            candidate.gain = splitGain.get(candidate.splitter);
            if candidate.gain > best.gain
              best = candidate;
            end
          end
        end
      end
    end
  end
end