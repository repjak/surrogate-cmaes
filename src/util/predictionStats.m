function s = predictionStats(y, ym, stat)
% s = predictionStats(y, ym, stat) count statistic for model prediction
%
% Input:
%   y    - original data
%   ym   - model prediction
%   stat - statistic to compute:
%            kendall  - Kendall' tau rank
%            mse      - mean square error
%            mzoe     - mean zero-one error
%            rankmse  - mean square error of ranks
%            rankmzoe - mean zero-one error of ranks
%            rde      - ranking difference error

  n = length(y);
  assert(length(ym) == n, 'Lengths of compared vectors has to be equal.')
  
  % compute ranking statistics
  if length(stat) > 3 && strcmp(stat(1:4), 'rank')
    [~, id] = sort(y);
    y = id(id);
    [~, id] = sort(ym);
    ym = id(id);
    stat = stat(5:end);
  end

  switch stat
    case 'mse'
      s = sum((ym - y).^2) / n;
    case 'mzoe'
      s = sum(ym ~= y) / n;
    case 'kendall'
      s = corr(ym, y, 'type', 'kendall');
    case 'rde'
      s = errRankMu(ym, y, floor(n/2));
    otherwise
      s = NaN;
  end

end