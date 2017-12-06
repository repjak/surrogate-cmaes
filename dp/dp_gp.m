% for reproducibility
rng('default');

phi = (1 + sqrt(5))/2;
xl = 10;
yl = xl/phi/2;
density = 100;
X = linspace(-xl, xl, density);
markerSize = 15;
fontSize = 30;
dp_path = '/home/jakub/Documents/vpdd2017active/img';

%% prior & posterior
exc = 1;
tries = 0;
while exc && tries < 10
  % prior
  covfunc = @covSEiso;
  ell = 2; sf = 1; hyp.cov = log([ell; sf]);
  meanfunc = @meanZero;
  likfunc = @likGauss; sn = 0.1; hyp.lik = log(sn);

  K = feval(covfunc, hyp.cov, X');
  yprior = mvnrnd(zeros(3, size(K, 1)), K);

  % posterior
  xs = [-7 -3 5 6];
  Ks = feval(covfunc, hyp.cov, xs');
  %ys = sin(xs);
  ys = mvnrnd(zeros(length(xs), 1), Ks);

  % axis manual;
  % xlim([-xl xl]);
  % ylim([-xl/phi/2 xl/phi/2]);

  hyp2 = minimize(hyp, @gp, -100, @infExact, [], covfunc, likfunc, xs', ys');
  K2 = feval(covfunc, hyp2.cov, X');
  [ymu, ys2, fmu, fs2, lp, post] = gp(hyp2, @infExact, meanfunc, covfunc, likfunc, xs', ys', X');
  %y = ymu + ys2.*randn(length(X), 1);

  KXXs = feval(covfunc, hyp2.cov, X', xs');
  KXsXs = feval(covfunc, hyp2.cov, xs');
  
  % Posterior covariance matrix
  kPost = K2 - KXXs/(KXsXs)*KXXs';

  try
    y = mvnrnd(repmat(ymu, 1, 5)', kPost);
    exc = 0;
  catch err
    exc = 1;
  end
  tries = tries + 1;
end

if exc
  close all;
  error('exception, exiting');
end

handle = figure();
set(handle,'PaperPositionMode','auto');

%ax1 = subplot(1, 3, 1);
plot(X, yprior(1,:), 'k-', X, yprior(2,:), 'k--', X, yprior(3,:), 'k-.');

% set(ax1, 'XLimMode', 'manual');
% set(ax1, 'YLimMode', 'manual');
% set(ax1, 'XLim', [-xl xl]);
% set(ax1, 'YLim', [-yl yl]);

xlim([-xl xl]);
ylim([-yl yl]);

set(gca, 'FontSize', fontSize);
set(gca, 'TickLabelInterpreter', 'latex');

xlabel('$x$', 'Interpreter', 'latex');
ylabel('$f(x)$', 'Interpreter', 'latex');
  
print(handle, '-depsc', fullfile(dp_path, 'gp_01'));

%ax2 = subplot(1, 3, 2);
plot(X, y(1,:), 'k-', X, y(2,:), 'k--', X, y(3,:), 'k-.', ...
  xs, ys, '+', 'Color', 'k', 'MarkerSize', markerSize);

% set(ax2, 'XLimMode', 'manual');
% set(ax2, 'YLimMode', 'manual');
% set(ax2, 'XLim', [-xl xl]);
% set(ax2, 'YLim', [-yl yl]);
xlim([-xl xl]);
ylim([-yl yl]);

set(gca, 'FontSize', fontSize);
set(gca, 'TickLabelInterpreter', 'latex');

xlabel('$x$', 'Interpreter', 'latex');
ylabel('$f(x)$', 'Interpreter', 'latex');

print(handle, '-depsc', fullfile(dp_path, 'gp_02'));


%% prediction with uncertainty

handle = figure();
set(handle,'PaperPositionMode','auto');

%ax3 = subplot(1, 3, 3);

% set(ax2, 'XLimMode', 'manual');
% set(ax2, 'YLimMode', 'manual');
% set(ax2, 'XLim', [-xl xl]);
% set(ax2, 'YLim', [-yl yl]);

plot(X, ymu);
hold on;

axis manual;
xlim([-xl xl]);
ylim([-yl yl]);

set(gca, 'FontSize', fontSize);
set(gca, 'TickLabelInterpreter', 'latex');

confplot(X, ymu, sqrt(ys2), sqrt(ys2), 'k');
hold on;
plot(xs, ys, '+', 'Color', 'k', 'MarkerSize', markerSize);

xlabel('$x$', 'Interpreter', 'latex');
ylabel('$f(x)$', 'Interpreter', 'latex');

% title('squared exponential (prediction)');

print(handle, '-depsc', fullfile(dp_path, 'gp_03'));

%print(handle, '-depsc', fullfile(dp_path, 'gp'));

close all;

%% probability of improvement

handle = figure();
set(handle, 'PaperPositionMode', 'auto');

confplot(X, ymu, sqrt(ys2), sqrt(ys2), 'k');
hold on;

plot(X, ymu, 'Color', 'k');
hold on;
plot(X, ymu - sqrt(ys2), 'Color', [0, 0.5, 0]);
plot(xs, ys, '+', 'Color', 'k', 'MarkerSize', markerSize);
xlim([-xl, xl]);
ylshift = 1;
ylim([-yl-ylshift, yl-ylshift]);

ymin = min(ymu);
ymin_ind = find(ymu == ymin, 1);
xmin = X(ymin_ind);
scale = 2.5;

x0 = xmin + 2;
[zmu, zs2] = gp(hyp2, @infExact, meanfunc, covfunc, likfunc, xs', ys', x0);

xpoi = linspace(-xl, xl, 1e4);
poi = normpdf(xpoi, 0, sqrt(zs2));
handle_poi = plot(xpoi, poi, 'Color', 'k');

T = (zmu - ymin)/scale;
xfill = xpoi(xpoi >= T);
handle_fill = fill([xfill, fliplr(xfill)], [poi(end-length(xfill)+1:end), ...
   repmat(0,1,length(xfill))], 'r');

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
xticklabels({'$x$'});
yticklabels({'$T$', '$\hat{f}(x)$'});
print(handle, '-depsc', fullfile(dp_path, 'gp_poi'));
close all;