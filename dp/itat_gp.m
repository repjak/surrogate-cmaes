% for reproducibility
rng('default');

phi = (1 + sqrt(5))/2;
xl = 10;
yl = xl/phi/2;
density = 100;
X = linspace(-xl, xl, density);
markerSize = 15;
fontSize = 25;
dp_path = '/home/jakub/Documents/dp/img/';

%% prior & posterior
exc = 1;
while exc
  handle = figure()%'Units', 'inches');
  set(handle,'PaperPositionMode','auto');
  %set(handle, 'PaperUnits','inches');
  %set(handle, 'Units','inches');
  %pos=get(handle,'Position');
  %set(handle, 'PaperSize', [6 5]);
  %set(handle, 'PaperPosition', [0 0 6 5]);

  axis manual;
  xlim([-xl xl]);
  ylim([-yl yl]);

  % prior
  covfunc = @covSEiso; ell = 2; sf = 1; hyp.cov = log([ell; sf]);
  meanfunc = @meanZero;
  likfunc = @likGauss; sn = 0.1; hyp.lik = log(sn);

  K = feval(covfunc, hyp.cov, X');
  %y = chol(K)'*gpml_randn(0.15, length(X), 1);% + exp(hyp.lik)*gpml_randn(0.2, length(X), 1);
  y = mvnrnd(zeros(3, size(K, 1)), K);

  plot(X, y(1,:), 'k-');
  hold on;
  plot(X, y(2,:), 'k--');
  plot(X, y(3,:), 'k-.');
  
  set(gca, 'FontSize', fontSize);
  xlabel('$x$', 'Interpreter', 'latex');
  ylabel('$f(x)$', 'Interpreter', 'latex');
  % title('squared exponential (prior)');

  print(handle, '-depsc', fullfile(dp_path, 'gp_01'));

  % posterior
  xs = [-7 -3 5 6];
  Ks = feval(covfunc, hyp.cov, xs');
  ys = mvnrnd(zeros(length(xs), 1), Ks);

  % axis manual;
  % xlim([-xl xl]);
  % ylim([-xl/phi/2 xl/phi/2]);

  %ell = 2; sf = 5; hyp.cov = log([ell; sf]);
  hyp2 = minimize(hyp, @gp, -100, @infExact, [], covfunc, likfunc, xs', ys');
  K2 = feval(covfunc, hyp2.cov, X');
  %post = infExact(hyp, {@meanZero}, covfunc, likfunc, xs', ys')
  [ymu ys2 fmu fs2 lp post] = gp(hyp2, @infExact, meanfunc, covfunc, likfunc, xs', ys', X');
  %y = ymu + ys2.*randn(length(X), 1);

  KXXs = feval(covfunc, hyp2.cov, X', xs');
  KXsXs = feval(covfunc, hyp2.cov, xs');
  sigma = exp(2*hyp.lik)*eye(max(size(xs))); 

  % Posterior covariance matrix
  kPost = K2 - KXXs/(KXsXs)*KXXs';

  try
    y = mvnrnd(repmat(ymu, 1, 5)', kPost);
    exc = 0;
  catch err
    exc = 1;
  end
end

handle = figure(); %'Units', 'inches');
%set(handle,'PaperPositionMode','manual');
%set(handle, 'PaperUnits','inches');
%set(handle, 'Units','inches');
set(handle,'PaperPositionMode','auto');


plot(X, y(1,:), 'k-');
hold on;
plot(X, y(2,:), 'k--');
plot(X, y(3,:), 'k-.');

%plot(X, ymu, 'Color', 'g');
plot(xs, ys, '+', 'Color', 'k', 'MarkerSize', markerSize);

set(gca, 'FontSize', fontSize);
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$f(x)$', 'Interpreter', 'latex');
% title('squared exponential (posterior)');

print(handle, '-depsc', fullfile(dp_path, 'gp_02'));


%% prediction with uncertainty

handle = figure();
set(handle,'PaperPositionMode','auto');

confplot(X, ymu, sqrt(ys2), sqrt(ys2), 'k');
hold on;
plot(xs, ys, '+', 'Color', 'k', 'MarkerSize', markerSize);

set(gca, 'FontSize', fontSize);
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$f(x)$', 'Interpreter', 'latex');
% title('squared exponential (prediction)');

print(handle, '-depsc', fullfile(dp_path, 'gp_03'));

close all;

%% probability of improvement

handle = figure();
set(handle, 'PaperPositionMode', 'auto');

plot(X, ymu, 'Color', 'k');
hold on;
plot(xs, ys, '+', 'Color', 'k', 'MarkerSize', markerSize);
xlim([-xl, xl]);
ylshift = 1;
ylim([-yl-ylshift, yl-ylshift]);

ymin = min(ymu);
ymin_ind = find(ymu == ymin, 1);
xmin = X(ymin_ind);
scale = 2.5;

x0 = xmin + 1.5;
[zmu, zs2] = gp(hyp2, @infExact, meanfunc, covfunc, likfunc, xs', ys', x0);

xpoi = linspace(-xl, xl, 1e4);
poi = normpdf(xpoi, 0, sqrt(zs2));
handle_poi = plot(xpoi, poi, 'Color', 'k');

T = (zmu - ymin)/scale;
xfill = xpoi(xpoi >= T);
handle_fill = fill([xfill, fliplr(xfill)], [poi(end-length(xfill)+1:end), ...
   repmat(0,1,length(xfill))], [0.8 0.8 0.8]);

plot([-xl, xl], [zmu zmu], 'Color', 'k');
plot([-xl, xl], [ymin ymin], 'Color', 'k');
plot([x0, x0], [(-yl-ylshift) (yl-ylshift)], 'Color', 'k');

t = hgtransform('Parent', gca);
s = hgtransform('Parent', gca);

set(handle_poi, 'Parent', t);
set(handle_fill, 'Parent', s);
T1 = makehgtform('zrotate', -pi/2, 'translate', [-zmu x0 0], 'scale', scale);
set(t, 'Matrix', T1);
set(s, 'Matrix', T1);

set(gca, 'FontSize', 20);
xlabel('');
ylabel('');
xticks([x0]);
yticks([ymin, zmu]);
set(gca, 'TickLabelInterpreter', 'latex');
xticklabels({'$x_\ast$'});
yticklabels({'$T$', '$f_\ast(x_\ast)$'});
print(handle, '-depsc', fullfile(dp_path, 'gp_poi'));