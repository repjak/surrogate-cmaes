function [freq, err] = modelStats(exptype, dims, funs)
%MODELSTATS Summary of this function goes here
%   Detailed explanation goes here

infodir = 'out/info/';
ninst = 15;

ndims = length(dims);
nfuns = length(funs);
freq = cell(nfuns, ndims);
err = cell(nfuns, ndims);

for dind = 1:length(dims)
  for find = 1:length(funs)
    for inst = 1:ninst
      suffix = sprintf('%d_%d_%d', dind, find, inst);
      gmodel_file = fullfile(infodir, [exptype '_gmodel.' suffix]);
      gorig_file = fullfile(infodir, [exptype '_gorig.' suffix]);
      err_file = fullfile(infodir, [exptype '_err.' suffix]);
      gmodel = dlmread(gmodel_file);
      gorig = dlmread(gorig_file);
      errs = dlmread(err_file);
      freq{find, dind} = [freq{find, dind} gorig/gmodel];
      err{find, dind} = [err{find, dind} mean(errs)];
    end
  end
end

end
