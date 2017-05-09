dp_path = '/home/jakub/Documents/dp/img/';

handle = figure('Units', 'inches');

set(handle,'PaperPositionMode','manual');
set(handle, 'PaperUnits','inches');
set(handle, 'Units','inches');
% pos=get(handle,'Position');
% set(handle, 'PaperSize', [5.5 5]);
% set(handle, 'PaperPosition', [0 0 5.5 5]);

fontSize = 15;

X = linspace(0, 30, 1e3);

ell=10;
sf=sqrt(1);
hyp.cov = log([ell; sf]);
K1 = feval(covfcn, hyp.cov, X', repmat(0, length(X), 1));

plot(X, K1(:,1), 'k--');

hold on;
ell=5;
sf=sqrt(1);
hyp.cov = log([ell; sf]);
K2 = feval(covfcn, hyp.cov, X', repmat(0, length(X), 1));

plot(X, K2(:,1), 'k-.');

hold on;
ell=5;
sf=sqrt(0.5);
hyp.cov = log([ell; sf]);
K3 = feval(covfcn, hyp.cov, X', repmat(0, length(X), 1));
plot(X, K3(:,1), 'k-');

xlabel('$x$', 'Interpreter', 'latex');
ylabel('$\mathrm{cov^{SE}}(0, x, \theta, l)$', 'Interpreter', 'latex');
set(gca, 'FontSize', fontSize);
legend({'$\theta=1.0, l=10$', '$\theta=1.0, l=5$', '$\theta=\frac{1}{2}, l=5$'}, 'Interpreter', 'latex');
axis square;

print(handle, '-depsc', fullfile(dp_path, 'covse_01'));
close all;