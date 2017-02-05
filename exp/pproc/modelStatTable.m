function modelStatTable(stats, varargin)
% modelStatTable(data, settings)
% Creates and prints table of model statistics.
%
% Input:
%   stats     - structure of model statistics
%   settings - pairs of property (string) and value or struct with 
%              properties as fields:
%
%     'DataNames'   - cell array of data names (e.g. names of algorithms)
%     'DataDims'    - dimensions of data
%     'DataFuns'    - functions of data
%     'Evaluations' - evaluations chosen to count
%     'Format'      - table format | ('tex', 'figure')
%     'Ranking'     - type of ranking (see help createRankingTable)
%                       'tolerant' - equal rank independence
%                       'precise'  - equal ranks shift following ranks
%                       'median'   - equal ranks replaced by medians of
%                                    shifted ranks (from 'precise')
%     'ResultFile'  - file containing resulting table
%     'Statistic'   - statistic of data | string or handle (@mean, @median)
%     'TableDims'   - dimensions chosen to count
%     'TableFuns'   - functions chosen to count
%
% Output:
%   rankTable - table of rankings
%   ranks     - rankings for each function and dimension
%
% See Also:
%   createRankingTable, speedUpPlot, speedUpPlotCompare, dataReady

  % initialization
  if nargin < 1 || isempty(stats) || ~isstruct(stats)
    help modelStatTable
    return
  end
  settings = settings2struct(varargin);

  % find statistics in stats
  statNames = fieldnames(stats);
  statValues = cellfun(@(x) stats.(x), statNames, 'UniformOutput', false);
  
  % parse settings
  [numOfModels, nFuns, nDims] = size(statValues{1});
  modelNames = defopts(settings, 'ModelNames', ...
    arrayfun(@(x) sprintf('model%02d', x), 1:numOfModels, 'UniformOutput', false));
  defaultDims = [2, 3, 5, 10, 20, 40];
  dataDims   = defopts(settings, 'DataDims', defaultDims(1:nDims));
  dataBBfuns = defopts(settings, 'DataFuns', 1:nFuns);
  tableFormat = defopts(settings, 'Format', 'disp');
  tableType   = defopts(settings, 'TableType', 'meanstd');
  dims    = defopts(settings, 'TableDims', dataDims);
  BBfunc  = defopts(settings, 'TableFuns', dataBBfuns);
  defResultFolder = fullfile('exp', 'pproc', 'tex');
  resultFolder = defopts(settings, 'ResultFolder', defResultFolder);
  
  % return only dimensions in the original data
  dimsId = ismember(dims, dataDims);
  funcId = ismember(BBfunc, dataBBfuns);
  dims = dims(dimsId);
  BBfunc = BBfunc(funcId);
  
  % prepare data
  switch tableType
    case 'meanstd'
      mark.mean = cellfun(@(x) strcmp(x(1:4), 'mean'), statNames);
      mark.std = cellfun(@(x) strcmp(x(1:3), 'std'), statNames);
      mainStatNames = cellfun(@(x) x(5:end), statNames(mark.mean), 'UniformOutput', false);
      % mark pairs of means and stds
      for s = 1:length(statNames)
        mark.stat(s) = find(~cellfun(@isempty, regexp(statNames{s}, mainStatNames)));
      end
      numOfTables = length(mainStatNames);
      % prepare table data
      for t = 1:numOfTables
        tableData(t).mean = stats.(['mean', mainStatNames{t}]);
        tableData(t).std  = stats.(['std',  mainStatNames{t}]);
      end
    otherwise
      warning('There is no table type: %s. Ending table printing.')
      help modelStatTable
      return
  end
  
  % print table
  switch tableFormat
    % display table
    case 'disp'
      for s = 1:length(statNames)
        fprintf('*** %s ***\n', statNames{s})
        dispTable(statValues{s}, BBfunc, dims, modelNames)
      end
      
    % prints table to latex file
    case {'tex', 'latex'}
      if ~exist(resultFolder, 'dir')
        mkdir(resultFolder)
      end
      for t = 1:numOfTables
        resultFile = fullfile(resultFolder, [mainStatNames{t}, '.tex']);
        FID = fopen(resultFile, 'w');
        % print table to tex file
        printTableTex(FID, tableData(t), mainStatNames{t}, BBfunc, dims, modelNames)
        fclose(FID);
        fprintf('Table written to %s\n', resultFile);
      end
  end

end

function dispTable(stat, func, dims, modelNames)
% display statistics in table
  for d = 1:length(dims)
    fprintf('  %dD  ', dims(d))
    fprintf('%s\n', strjoin(modelNames))
    for f = 1:length(func)
      fprintf(' f%02d ', func(f))
      for m = 1:size(stat, 1)
        if isnan(stat(m, f, d))
          fprintf('    ---  ')
        elseif stat(m, f, d) < 0
          fprintf('  %0.3f ', stat(m, f, d))
        else
          fprintf('%s%0.3f ', ones(1, 3-max(0, floor(log10(stat(m, f, d)))))*' ', stat(m, f, d))
        end
      end
      fprintf('\n')
    end
    fprintf('\n')
  end
end

function printTableTex(FID, tableData, mainStatName, func, dims, modelNames)
% Prints table to file FID

  nFunc = length(func);
  nDims = length(dims);
  nModels = length(modelNames);
  
  % specifications of statistic
  [statName, boldFunc] = statSpecification(mainStatName);
  
  % table headers
  fprintf(FID, '\\begin{table}\n');
  fprintf(FID, '\tiny\n');
  fprintf(FID, '\\centering\n');
  fprintf(FID, '\\resizebox*{\\textwidth}{\\textheight}{');
  fprintf(FID, '\\begin{tabular}[pos]{ l %s }\n', repmat(' | c', 1, nModels));
  fprintf(FID, '\\hline\n');
  % dimension section
  for d = 1:nDims
    fprintf(FID, ' \\textbf{%dD} & %s \\\\\n', dims(d), strjoin(modelNames, ' & '));
    fprintf(FID, '\\hline\n');
    % print function row
    for f = 1:nFunc
      fprintf(FID, ' f%d', func(f));
      fprintf(FID, '%s ', boldBestMeanStd(tableData.mean(:, f, d), tableData.std(:, f, d), boldFunc));
      fprintf(FID, '\\\\\n');
    end
    fprintf(FID, '\\hline\n');
  end
  fprintf(FID, '\\hline\n');
  fprintf(FID, '\\end{tabular}\n');
  fprintf(FID, '}\n');
  % dimension numbers 
  dimString = num2str(dims, ', %d');
  dimString = dimString(2:end);
  % caption printing
  fprintf(FID, '\\vspace{1mm}\n');
  fprintf(FID, ['\\caption{Means and standard deviations of %s of offline testing on 10 chosen generations ', ...
                'from %d benchmark functions and dimensions D = \\{%s\\}. ', ...
                'The best achieved values in each function are given in bold.}\n'], ...
                statName, nFunc, dimString);
               
  fprintf(FID, '\\label{tab:%s}\n', mainStatName);
  fprintf(FID, '\\end{table}\n');
  
end

function [statisticName, boldFunc] = statSpecification(statName)
% [statisticName, boldFunc] = statSpecification(statName) returns specific
% values for each statistic
  switch statName
    % ranking difference error
    case 'rde'
      statisticName = 'RDE';
      boldFunc = @min;
    % Kendall's tau rank
    case 'kendall'
      statisticName = 'Kendall''s tau rank';
      boldFunc = @max;
    % mean zero-one error of rank (Hamilton's distance)
    case 'rankmzoe'
      statisticName = 'MZOE';
      boldFunc = @min;
    otherwise
      statisticName = upper(statName);
      boldFunc = @min;
  end
end

function str = boldBestMeanStd(dataMean, dataStd, stat)
% returns row of MEAN ± STD, where the statistical best value is bolded
  nData = length(dataMean);
  bestId = (dataMean == stat(dataMean));
  msVals = cell(1, nData);
  for d = 1:nData
    if bestId(d)
      msVals{d} = sprintf('$ \\mathbf{%0.2f \\pm %0.2f} $', dataMean(d), dataStd(d));
    else
      msVals{d} = sprintf('$ %0.2f \\pm %0.2f $', dataMean(d), dataStd(d));
    end
  end
  str = [' & ', strjoin(msVals, ' & ')];
end
