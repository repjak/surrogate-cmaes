[X, Y] = meshgrid(linspace(-2, 2, 1e2), linspace(-2, 2, 1e2));
opt = [10 10];
q = @(x) norm(x-opt)^2;
sigma = 1;
lambda = 100;
mu = ceil(0.30 * lambda);
w = [1/mu 1/mu repmat((1-1/mu-1/mu)/(mu-2), 1, mu-2)];
C = diag([sigma sigma]);

handle = figure('Units', 'inches');

dp_path = '/home/jakub/Documents/dp/img/';

% objective contours
contour(X, Y, reshape(sqrt((X-opt(1)).^2+(Y-opt(2)).^2), length(X), length(Y)), 'LineColor', [0.75 0.75 0.75], 'LevelStep', 0.75);
hold on;
axis equal;

% CM contour
contour(X, Y, reshape(mvnpdf([X(:) Y(:)], 0, C), length(X), length(Y)), repmat(mvnpdf(sqrt([0 sigma]), 0, C), 1, 2), 'LineColor', 'k');
%plot([0], sqrt([sigma]), 'o', 'Color', 'r')

% mean
plot([0], [0], '+', 'Color', 'black', 'MarkerSize', 6, 'MarkerFaceColor', 'k');

set(gca, 'XTick', []);
set(gca, 'YTick', []);
set(gca, 'XLabel', []);
set(gca, 'YLabel', []);

set(handle,'PaperPositionMode','manual');
set(handle, 'PaperUnits','inches');
set(handle, 'Units','inches');
pos=get(handle,'Position');
set(handle, 'PaperSize', [3 3]);
set(handle, 'PaperPosition', [-0.41 -0.41 3.7 3.7]);

print(handle, '-dpdf', '-r1500', fullfile(dp_path, 'cma_ellipse_01'));



clf;

% objective contours
h = contour(X, Y, reshape(sqrt((X-opt(1)).^2+(Y-opt(2)).^2), length(X), length(Y)), 'LineColor', [0.75 0.75 0.75], 'LevelStep', 0.75);
hold on;
axis equal;

% mean
plot([0], [0], '+', 'Color', 'black', 'MarkerSize', 6, 'MarkerFaceColor', 'k');

B = eye(2, 2);
D = [0.25 0; 0 1.5^2];
Cn = B*D*B';

% CM contour
contour(X, Y, reshape(mvnpdf([X(:) Y(:)], 0, Cn), length(X), length(Y)), repmat(mvnpdf(sqrt(D(1,1)) * B(:,1), 0, Cn), 1, 2), 'LineColor', 'k', 'ShowText', 'Off');

set(gca, 'XTick', []);
set(gca, 'YTick', []);
set(gca, 'XLabel', []);
set(gca, 'YLabel', []);

print(handle, '-dpdf', '-r1500', fullfile(dp_path, 'cma_ellipse_02'));



clf;
% objective contours
contour(X, Y, reshape(sqrt((X-opt(1)).^2+(Y-opt(2)).^2), length(X), length(Y)), 'LineColor', [0.75 0.75 0.75], 'LevelStep', 0.75);
hold on;
axis equal;

angle = pi/6;
B = [cos(angle) sin(angle); -sin(angle) cos(angle)];
sigma = 5
D = [0.25 0; 0 1.5^2];
Cn = B*D*B';

% mean
plot([0], [0], '+', 'Color', 'black', 'MarkerSize', 6, 'MarkerFaceColor', 'k');

% CM contour
contour(X, Y, reshape(mvnpdf([X(:) Y(:)], 0, Cn), length(X), length(Y)), repmat(mvnpdf(sqrt(D(1,1)) * B(:,1), 0, Cn), 1, 2), 'LineColor', 'k', 'ShowText', 'Off');


set(gca, 'XTick', []);
set(gca, 'YTick', []);
set(gca, 'XLabel', []);
set(gca, 'YLabel', []);

print(handle, '-dpdf', '-r1500', fullfile(dp_path, 'cma_ellipse_03'));