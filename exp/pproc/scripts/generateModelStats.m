dims = [10];
funs = 1:24;
datanames = {'ADA-KL', 'ADA-Kendall', 'ADA-RD'};

colors = [[100, 149, 237]; ...
  [255, 0, 0]; ...
  [154, 205, 50]];

if (exist('modelstats.mat', 'file'))
  load('modelstats.mat');
else 
  freqs = cell(1,3);
  errs = cell(1,3);
  
  i = 1;
  for exptype = {'17_kl', '18_kendall', '19_rankdiff'};
    [freq, err] = modelStats(exptype{:}, dims, funs);
    freqs{i} = freq;
    errs{i} = err;
    
    i = i + 1;
  end
  save('modelstats.mat', 'errs', 'freqs');
end

for dim = dims
  printModelStats(datanames, funs, dim, freqs, errs, colors);
  close all;
end
