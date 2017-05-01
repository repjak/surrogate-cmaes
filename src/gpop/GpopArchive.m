classdef GpopArchive < Archive
  methods (Access = private, Static)
    function d = sq_dist(X, Y)
      % compute euclidean distance of row vectors in @X and @Y
      if (size(X, 1) == 1)
        d = [];
        for y = Y'
          d = [d; sqrt(sum((X - y').^2, 2))];
        end
      else
        assert(isequal(size(X), size(Y)), 'Dimensions don''t match');
        d = sqrt(sum((X - Y).^2, 2));
      end
    end
  end

  methods (Access = public)
    function obj = GpopArchive(dimension)
      obj@Archive(dimension);
    end

    function idx = isMember(obj, X, varargin)
      % test membership of a point in the archived data, optionally within a tolerated distance
      % given in an extra argument
      if (~isempty(obj.X))
        if (nargin > 2)
          idx = [];
          for x = X'
            idx = [idx; sum(obj.sq_dist(x', obj.X) < varargin{1}) > 0];
          end
        else
          idx = ismember(X, obj.X, 'rows');
        end
      else
        idx = true(size(X,1),1);
      end % if
    end % function

    function [X, y] = getNearData(obj, n, x)
      % return @n closest points to point @x in euclidean metrics
      nData = length(obj.y);

      if (nData <= n)
        X = obj.X;
        y = obj.y;
        return;
      end

      % euclidean distance from x
      dist = obj.sq_dist(obj.X, repmat(x, nData, 1));

      % sort all data by distance from x
      [sorted, sortedIdx] = sort(dist);
      idx = sortedIdx(1:n, :);

      X = obj.X(idx, :);
      y = obj.y(idx);
    end

    function [X, y, idx] = getRecentData(obj, n)
      % return up to @n most recently added points
      nData = length(obj.y);
      X = []; y = []; idx = [];

      if (nData == 0)
        return;
      else
        X = obj.X(max(1, nData - n + 1):nData, :);
        y = obj.y(max(1, nData - n + 1):nData, :);
        return;
      end
    end
  end % methods

end % classdef
