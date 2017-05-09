xl = 10;
yl = [-5 5];
yrange = yl(2) - yl(1);
density = 1e4;
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

axis manual;
xlim([-xl xl]);
ylim([-xl/phi/2 xl/phi/2]);


%% making prediction

covfunc = {@covSum, {@covSEiso, @covNoise}};
%covfunc = {@covSEiso};
ell = 1; sf = 2; snoise = 1; hyp.cov = [log([ell; sf; snoise])];
meanfunc = @meanZero;
likfunc = @likGauss; sn = 10; hyp.lik = log(sn);

xs = [-6:0.5:-1 3:0.5:5];
Ks = feval(covfunc{:}, hyp.cov, xs');
ys = mvnrnd(zeros(1, size(Ks, 1)), Ks);
%ys = ys + 0.2*randn(1, length(xs));

covfunc = {@covSEiso};
ell = 1; sf = 0.5;
hyp.cov = [log([ell; sf])];
hyp2 = minimize(hyp, @gp, -100, @infExact, [], covfunc, likfunc, xs', ys');
%K2 = feval(covfunc{:}, hyp2.cov, X');
[ymu ys2] = gp(hyp2, @infExact, meanfunc, covfunc, likfunc, xs', ys', X');

fill([X, fliplr(X)], [ymu' + sqrt(ys2') fliplr(ymu' - sqrt(ys2'))], ...
  0.9*[1 1 1], 'EdgeColor', [1 1 1]);
hold on;

%confplot(X, ymu, sqrt(ys2));
plot(X, ymu, 'k-');

hold on;
plot(xs, ys, 'r.', 'MarkerSize', markerSize);

x0 = -7;
[z zmu] = gp(hyp2, @infExact, meanfunc, covfunc, likfunc, xs', ys', x0)

% problem.solver = 'fminunc';
% problem.objective = @(x)gp(hyp2, @infExact, meanfunc, covfunc, likfunc, xs', ys', [x]);
% problem.x0 = [-xl,xl];
% problem.options = optimoptions(@fminunc,'Algorithm','quasi-newton','Display','off');

ymin = min(ymu);
ymin_ind = find(ymu == ymin, 1);
xmin = X(ymin_ind);

zdom = linspace(-xl, xl, 1e3);
zpred = yrange*normpdf(zdom, 0, 0.8*sqrt(zmu));
hpdf = plot(zdom, zpred, 'Color', 'k');

T = (z-ymin);
xfill = zdom(zdom >= T);
hfill = fill([xfill, fliplr(xfill)], [zpred(end-length(xfill)+1:end), repmat([0],1,length(xfill))], 'g');

t = hgtransform('Parent', gca);
s = hgtransform('Parent', gca);

set(hpdf, 'Parent', t);
set(hfill, 'Parent', s);
T1 = makehgtform('zrotate', -pi/2, 'translate',[-z x0 0]);
set(t, 'Matrix', T1);
set(s, 'Matrix', T1);

plot(x0, z, 'k.', 'MarkerSize', markerSize, 'Color', 'k');
text(-xl, z, '$\hat{f}_{x}$', 'Interpreter', 'latex', 'FontSize', fontSize, 'HorizontalAlignment', 'right');
text(x0, yl(1), '$x$', 'Interpreter', 'latex', 'FontSize', fontSize, 'VerticalAlignment', 'top');
text(-xl, ymin, '$T$', 'Interpreter', 'latex', 'FontSize', fontSize, 'HorizontalAlignment', 'right');
plot([-xl, x0+yrange*normpdf(0, 0, 0.8*sqrt(zmu))], [z z], 'Color', 'k');
plot([-xl, xmin], [ymin ymin], 'Color', 'k');
line([x0 x0], yl, 'Color', 'k');

line([-xl -xl], yl, 'Color', 'k');
line([xl xl], yl, 'Color', 'k');


ylim(yl);

set(gca, 'xtick', []);
set(gca, 'ytick', []);
xlabel('');
ylabel('');
print(handle, '-dpdf', '-r1500', fullfile(dp_path, 'gp_poi'));