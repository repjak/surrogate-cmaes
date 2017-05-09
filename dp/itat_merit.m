xl = 10;
%yl = [-5 5];
%yrange = yl(2) - yl(1);
density = 1000;
X = linspace(-xl, xl, density);
phi = (1 + sqrt(5))/2;
markerSize = 20;
fontSize = 20;
dp_path = '/home/jakub/Documents/itat_slides/img/';

handle = figure('Units', 'inches');

set(handle,'PaperPositionMode','manual');
set(handle, 'PaperUnits','inches');
set(handle, 'Units','inches');
pos=get(handle,'Position');
set(handle, 'PaperSize', [6 5]);
set(handle, 'PaperPosition', [0 0 6 5]);

% axis manual;
% xlim([-xl xl]);
% ylim([-xl/phi/2 xl/phi/2]);


alpha = [0 1 4];

ftarget = @(x) sin(x);

Xtrain = [-6 0 4 6 6.1];
Ytrain = ftarget(Xtrain);

covfunc = {@covSum {@covSEiso, @covNoise}};
ell = 1; sf = 1; snoise = 1;
hyp.cov = [log([ell; sf; snoise])];
likfunc = @likGauss; sn = 10; hyp.lik = log(sn);
meanfunc = @meanZero;
hyp2 = minimize(hyp, @gp, -100, @infExact, [], covfunc, likfunc, Xtrain', Ytrain');
[ymu ys2] = gp(hyp2, @infExact, meanfunc, covfunc, likfunc, Xtrain', Ytrain', X');

fill([X fliplr(X)], [ymu' + ys2' fliplr(ymu' - ys2')], ...
  0.8*[1 1 1], 'EdgeColor', 0.8*[1 1 1]);
hold on;

plot(X, [ftarget(X); ymu'; ymu' - 4*ys2']);

plot(Xtrain, Ytrain, 'r.', 'MarkerSize', markerSize);

% x0 = -7;
% [z zmu] = gp(hyp2, @infExact, meanfunc, covfunc, likfunc, xs', ys', x0)
% 
% zdom = linspace(-xl, xl, 1e3);
% zpred = yrange*normpdf(zdom, 0, 0.8*sqrt(zmu));
% hpdf = plot(zdom, zpred, 'Color', 'k');
% 
% t = hgtransform('Parent', gca);
% set(hpdf, 'Parent', t);
% T1 = makehgtform('zrotate', -pi/2, 'translate',[-z x0 0]);
% set(t, 'Matrix', T1);
% plot(x0, z, 'k.', 'MarkerSize', markerSize, 'Color', 'k');
% t_h1 = text(-xl, z, '$\hat{y}_{N+1}$', 'Interpreter', 'latex', 'FontSize', fontSize, 'HorizontalAlignment', 'right');
% t_h2 = text(x0, yl(1), '$x_{N+1}$', 'Interpreter', 'latex', 'FontSize', fontSize, 'VerticalAlignment', 'top');
% %plot([x0, x0+10*normpdf(0, 0, sqrt(zmu))], [z z], 'Color', 'k');
% plot([-xl, x0+yrange*normpdf(0, 0, 0.8*sqrt(zmu))], [z z], 'Color', 'k');
% line([x0 x0], yl, 'Color', 'k');
% 
% ylim(yl);
% 
% set(gca, 'xtick', []);
% set(gca, 'ytick', []);
% xlabel('');
% ylabel('');
% print(handle, '-dpdf', '-r1500', fullfile(dp_path, 'gp_merit_01'));