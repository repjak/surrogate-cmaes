classdef TreeModel < AbstractModel
  properties (Constant, Access = private)
    nodeTemplate = struct(... % template for nodes
        'parent', 0, ...
        'left', 0, ...
        'right', 0, ...
        'splitter', [], ...
        'predictor', [], ...
        'X', [], ...
        'y', []);
  end
  
  properties %(Access = protected)
    nodes % list of nodes
    nNodes % nodes count
    minGain % minimum gain to split
    minLeafSize % minimum number of data examples in leaf
    minParentSize % minimum number of data examples in parent
    maxDepth % maximum depth of the tree
    splits % generators for split functions
    splitGain % evaluator for split functions
    predictorFunc % function which creates a model in leaf
    pruning % grows a full tree then prunes, otherwise prunes during splitting
    objectiveFunc
  end
  
  methods
    function obj = TreeModel(modelOptions, xMean)
      % constructor
      obj = obj@AbstractModel(modelOptions, xMean);
      
      % specific model options
      obj.minGain = defopts(modelOptions, 'minGain', 1e-3);
      obj.minLeafSize = defopts(modelOptions, 'minLeafSize', 2);
      obj.minParentSize = defopts(modelOptions, 'minParentSize', 10);
      obj.maxDepth = defopts(modelOptions, 'maxDepth', inf);
      obj.pruning = defopts(modelOptions, 'pruning', false);
      obj.objectiveFunc = defopts(modelOptions, 'objectiveFunc', @immse);
      obj.predictorFunc = defopts(modelOptions, 'predictorFunc', ...
        @(xMean) ConstantModel(struct, xMean));
      obj.splits = defopts(modelOptions, 'splits', ...
        {AxisSplit(1, 1, false)});
      obj.splitGain = defopts(modelOptions, 'splitGain', ...
        GenericModelTreesplitGain(@gainMse, obj.predictorFunc));
    end
    
    function nData = getNTrainData(obj)
      % returns the required number of data for training the model
      nData = obj.minLeafSize;
    end

    function obj = trainModel(obj, X, y, xMean, generation)
      % train the model based on the data (X,y)
      obj.trainGeneration = generation;
      obj.trainMean = xMean;
      obj.dataset.X = X;
      obj.dataset.y = y;
      
      obj.nNodes = 0;
      initialSize = min(2^obj.maxDepth, 2 * round(size(X,1) / obj.minLeafSize));
      obj.nodes = repmat(TreeModel.nodeTemplate, initialSize, 1);
      iNodeRoot = obj.addNode();
      obj.trainModelRecursive(X, y, iNodeRoot, 0);
    end

    function [y, sd2] = modelPredict(obj, X)
      % predicts the function values in new points X
      [y, sd2] = obj.modelPredictRecursive(X, 1);
    end
    
    function prune(obj, X, y)
      obj.pruneRecursive(X, y, 1);
    end
  end
  
  methods (Access = private)
    function trainModelRecursive(obj, X, y, iNode, depth)
      [N, D] = size(X);
      if depth < obj.maxDepth && N >= obj.minParentSize && N >= 2 * obj.minLeafSize && length(unique(y)) >= 2
        best = struct('gain', -inf);
        obj.splitGain.reset(X, y);
        for iSplit = 1:size(obj.splits, 2)
          split = obj.splits{iSplit};
          split.reset(X, y);
          candidate = split.get();
          idx = candidate.splitter(X);
          if sum(idx) >= obj.minLeafSize && sum(~idx) >= obj.minLeafSize
            if candidate.gain > best.gain
              best = candidate;
            end
          end
        end
        if best.gain >= obj.minGain
          obj.nodes(iNode).splitter = best.splitter;
          idx = best.splitter(X);
          
          left = struct('idx', idx, 'iNode', obj.addNode());
          obj.nodes(iNode).left = left.iNode;
          obj.nodes(left.iNode).parent = iNode;
          obj.trainModelRecursive(X(left.idx, :), y(left.idx, :), left.iNode, depth+1);
          
          right = struct('idx', ~idx, 'iNode', obj.addNode());
          obj.nodes(iNode).right = right.iNode;
          obj.nodes(right.iNode).parent = iNode;
          obj.trainModelRecursive(X(right.idx, :), y(right.idx, :), right.iNode, depth+1);
          
          return;
        end
      end
      if isempty(obj.nodes(iNode).predictor)
        obj.nodes(iNode).predictor = obj.predictorFunc(mean(X));
        obj.nodes(iNode).predictor.trainModel(X, y, mean(X), 0);
        obj.nodes(iNode).X = X;
        obj.nodes(iNode).y = y;
      end
    end
    
    function [y, sd2] = modelPredictRecursive(obj, X, iNode)
      if isempty(obj.nodes(iNode).splitter)
        [y, sd2] = obj.nodes(iNode).predictor.modelPredict(X);
      else
        y = zeros(size(X, 1), 1);
        sd2 = zeros(size(X, 1), 1);
        
        left = struct('idx', obj.nodes(iNode).splitter(X), 'iNode', obj.nodes(iNode).left);
        if any(left.idx)
          [y(left.idx), sd2(left.idx)] = obj.modelPredictRecursive(X(left.idx, :), left.iNode);
        end
        
        right = struct('idx', ~left.idx, 'iNode', obj.nodes(iNode).right);
        if any(right.idx)
          [y(right.idx), sd2(right.idx)] = obj.modelPredictRecursive(X(right.idx, :), right.iNode);
        end
      end
    end
    
    function [y, sd2] = modelPredictIterative(obj, X)
      y = zeros(size(X, 1), 1);
      sd2 = zeros(size(X, 1), 1);
      for i = 1:size(X, 1)
        iNode = 1;
        while ~isempty(obj.nodes(iNode).splitter)
          if obj.nodes(iNode).splitter(X(i, :))
            iNode = obj.nodes(iNode).left;
          else
            iNode = obj.nodes(iNode).right;
          end
        end
        [y(i), sd2(i)] = obj.nodes(iNode).predictor.modelPredict(X(i, :));
      end
    end
    
    function pruneRecursive(obj, X, y, iNode)
      if obj.nodes(iNode).left > 0 && obj.nodes(iNode).right > 0 ...
        && obj.nodes(obj.nodes(iNode).left) == 0 ...
        && obj.nodes(obj.nodes(iNode).right) == 0
        % internal node having two leaves
        % consider making this node a leaf
        yPred = zeros(size(X, 1), 1);
        
        left = struct('idx', obj.nodes(iNode).splitter(X), 'iNode', obj.nodes(iNode).left);
        if any(left.idx)
          [yPred(left.idx)] = obj.modelPredictRecursive(X(left.idx, :), left.iNode);
        end
        
        right = struct('idx', ~left.idx, 'iNode', obj.nodes(iNode).right);
        if any(right.idx)
          [yPred(right.idx)] = obj.modelPredictRecursive(X(right.idx, :), right.iNode);
        end
        
        objective = obj.objectiveFunc(y, yPred);
        
        current = struct;
        current.X = [obj.nodes(left.iNode).X; obj.nodes(right.iNode).X];
        current.y = [obj.nodes(left.iNode).y; obj.nodes(right.iNode).y];
        current.xMean = mean(X);
        predictorNew = obj.predictorFunc(current.xMean);
        predictorNew.trainModel(current.X, current.y, current.xMean, 0);
        objectiveNew = obj.objectiveFunc(y, predictorNew.modelPredict(X));
        
        if objectiveNew <= objective
          % remove children
          obj.nodes(left.iNode) = TreeModel.nodeTemplate;
          obj.nodes(right.iNode) = TreeModel.nodeTemplate;
          obj.nodes(iNode).left = 0;
          obj.nodes(iNode).right = 0;
          obj.nodes(iNode).predictor = predictorNew;
          obj.nodes(iNode).X = X;
          obj.nodes(iNode).y = y;
        end
      end
      if obj.nodes(iNode).left == 0 && obj.nodes(iNode).right == 0
        % try to replace complex model in leaf with constant model
        objective = obj.objectiveFunc(y, obj.nodes(iNode).predictor.modelPredict(X));
        
        xMean = mean(obj.nodes(iNode).X);
        predictorNew = ConstantModel(struct, xMean);
        predictorNew.trainModel(obj.nodes(iNode).X, obj.nodes(iNode).y, xMean, 0);
        objectiveNew = obj.objectiveFunc(y, predictorNew.modelPredict(X));
        
        if objectiveNew <= objective
          obj.nodes(iNode).predictor = predictorNew;
        end
      end
    end
    
    function iNode = addNode(obj)
      obj.nNodes = obj.nNodes + 1;
      if size(obj.nodes, 1) < obj.nNodes
        obj.children(2 * obj.nNodes, 1) = TreeModel.nodeTemplate;
      end
      iNode = obj.nNodes;
    end
  end
  
end