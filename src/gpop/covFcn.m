function Kmn = covFcn(Xm, Xn, theta)
% A covariance function defined as sum of SEard, constant
% shift and a white noise:
%
% Kmn(x_i, x_j) = sf2 * exp(-(x_i - x_q)'*diag(1./ell)*(x_i - x_j)/2)
%                 + sc2
%                 + sn2*delta(x_i, x_j)
%
% where ell = (ell_1^2,...ell_D^2) is vector of ARD parameters, D is
% dimension of input space, sf2 is signal variance, sc2 is constant shift
% and sn2 is added white noise.
%
% Theta is a vector of positive numbers:
% theta = [ log(ell_1)
%           log(ell_2)
%            .
%           log(ell_D)
%           log(sqrt(sf2))
%           log(sqrt(sc2))
%           log(sqrt(sn2)) ]
%
  assert(size(Xm, 2) == size(Xn, 2), 'Dimensions must agree');
  dim = size(Xm, 2);
  
  % tolerance on squared distance of vectors when their considered equal
  tol = 1e-9;

  assert(all(size(theta) == [dim+3 1]), ['theta must be a row vector of ''' dim + 3 ''' hyperparameters']);

  ell = exp(theta(1:dim));
  sf2 = exp(2*theta(dim+1));
  sc2 = exp(2*theta(dim+2));
  sn2 = exp(2*theta(dim+3));

  Kmn = sq_dist(diag(1./ell) * Xm', diag(1./ell) * Xn');
  Kmn = sf2 * exp(-Kmn/2);
  Kmn = bsxfun(@plus, Kmn, sc2);
  
  delta = bsxfun(@lt, sq_dist(Xm', Xn'), tol^2);
  Kmn(delta) = Kmn(delta) + sn2;
end


