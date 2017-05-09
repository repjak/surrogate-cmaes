% for reproducibility
rng('default');

xl = 5;
density = 15;
[X1, X2] = meshgrid(linspace(-xl, xl, density));
[U1, U2] = meshgrid(linspace(-xl, xl, density*3));
X = [X1(:) X2(:)];
markerSize = 15;
fontSize = 20;
phi = (1 + sqrt(5))/2;
dp_path = '/home/jakub/Documents/dp/img/';

handle = figure('Units', 'inches');
colormap gray;
%set(handle,'PaperPositionMode','manual');
%set(handle, 'PaperUnits','inches');
% set(handle, 'Units','inches');
% pos=get(handle,'Position');
%set(handle, 'PaperSize', [6 5]);
%set(handle, 'PaperPosition', [0 0 6 5]);


axis manual;
xlim([-xl xl]);
ylim([-xl xl]);
zlim([-xl xl]);


covfunc = @sqexpard;
ell = [1; 1];
sf = 1;
hyp.cov = log([sf; ell]);

%K = feval(covfunc, hyp.cov, X);
K = feval(covfunc, X, X, hyp.cov);
L = chol(K);
y = L'*gpml_randn(0.15, size(X, 1), 1);

v = griddata(X1, X2, reshape(y, density, density), U1, U2, 'cubic');
mesh(U1, U2, v);
set(gca, 'FontSize', fontSize);
set(gca, 'TickLabelInterpreter', 'latex');

xlabel('$x$', 'Interpreter', 'latex');
ylabel('$y$', 'Interpreter', 'latex');
print(handle, '-depsc', fullfile(dp_path, 'gp_ard_01'));


handle = figure('Units', 'inches');
colormap gray;
%set(handle,'PaperPositionMode','manual');
%set(handle, 'PaperUnits','inches');
% set(handle, 'Units','inches');
% pos=get(handle,'Position');
%set(handle, 'PaperSize', [6 5]);
%set(handle, 'PaperPosition', [0 0 6 5]);

ell = [2; 0.5];
sf = 1;
hyp.cov = log([sf; ell]);

%K = feval(covfunc, hyp.cov, X);
K = feval(covfunc, X, X, hyp.cov);
L = chol(K);
y = L'*gpml_randn(0.15, size(X, 1), 1);

v = griddata(X1, X2, reshape(y, density, density), U1, U2, 'cubic');

hold off;
mesh(U1, U2, v);
set(gca, 'FontSize', fontSize);
set(gca, 'TickLabelInterpreter', 'latex');

xlabel('$x$', 'Interpreter', 'latex');
ylabel('$y$', 'Interpreter', 'latex');
print(handle, '-depsc', fullfile(dp_path, 'gp_ard_02'));
%close all;