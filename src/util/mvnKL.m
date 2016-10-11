function divKL = mvnKL(m1, C1, m2, C2)
  % Compute Kullback-Leibler divergence of two multivariate
  % distributions given by means 'm1' and 'm2', respectively
  % and positive semidefinite covariance matrices 'C1' and 'C2',
  % respectively.
  % Return NaN if any of the matrices is not positive
  % definite.

  assert(length(m1) == length(m2) && all(size(C1) == size(C2)), 'Dimensions don''t match');

  [L1, p1] = chol(C1, 'lower');
  if p1 > 0, return NaN, end

  logdetC1 = GenerationsUpdaterKL.cov_logdet(L1);

  [L2, p2] = chol(C2, 'lower');
  if p2 > 0, return NaN, end

  logdetC2 = GenerationsUpdaterKL.cov_logdet(L2);
  invC2 = GenerationsUpdaterKL.cov_inv(L2);

  dim = length(m1);
  trprod = invC2(:)'*C1(:);
  mdiff = m2 - m1;

  divKL = 0.5 * (trprod + mdiff' * invC2 * mdiff - dim + (logdetC2 - logdetC1));
  divKL = max(0, divKL);
end

function cov_logdet = cov_logdet(L)
  % compute log(det(C)) for a positive semidefinite matrix given by
  % a lowertriangular matrix from a LU factorization of C
  cov_logdet = 2 * sum(log(diag(L)));
end

function cov_inv = cov_inv(L)
  % compute covariance matrix inverse from its lower triangular LU-factor
  invL = L\eye(size(L, 1));
  cov_inv = invL' * invL;
end
