[X, Y] = meshgrid(linspace(-3, 3, 1e2), linspace(-3, 3, 1e2));
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

pop = mvnrnd([0, 0], diag([sigma^2, sigma^2]), lambda); plot(pop(:,1), pop(:,2), '.', 'Color', 'b', 'MarkerFaceColor', 'b', 'MarkerSize', 4.0)
% population
plot(pop(:,1), pop(:,2), 'o', 'Color', 'black', 'MarkerFaceColor', 'k', 'MarkerSize', 1.5)

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

print(handle, '-dpdf', '-r1500', fullfile(dp_path, 'cma_est_a_01'));

f = [];
for x = pop'
  f(end+1) = q(x');
end



clf;

% objective contours
h = contour(X, Y, reshape(sqrt((X-opt(1)).^2+(Y-opt(2)).^2), length(X), length(Y)), 'LineColor', [0.75 0.75 0.75], 'LevelStep', 0.75);
hold on;
axis equal;

% mean
plot([0], [0], '+', 'Color', 'black', 'MarkerSize', 6, 'MarkerFaceColor', 'k');

[s, i] = sort(f);
pop = pop(i(1:mu),:);

%plot(pop(:,1), pop(:,2), 'o', 'Color', 'black', 'MarkerFaceColor', 'g', 'MarkerSize', 1.5)

for x = pop'
  y = [0 x(1)];
  z = [0 x(2)];
  %quiver(z(1), y(1), z(2)-z(1), y(2)-y(1), 0, 'Color', 'k');
  arrow(y, z);
  hold on;
end

% update Cn
Cn = zeros(2, 2);
for i = 1:mu
  Cn = Cn + w(i)*pop(i,:)'*pop(i,:);
end
%Cn = (1/(sigma^2))*Cn;

[B, D] = eig(Cn);
b1 = B(:,1);
b2 = B(:,2);
d = sqrt([D(1,1) D(2,2)]);

% CM contour
contour(X, Y, reshape(mvnpdf([X(:) Y(:)], 0, C), length(X), length(Y)), repmat(mvnpdf(sqrt([0 sigma]), 0, C), 1, 2), 'LineStyle', '--', 'LineColor', 'k', 'ShowText', 'Off');
contour(X, Y, reshape(mvnpdf([X(:) Y(:)], 0, Cn), length(X), length(Y)), repmat(mvnpdf(d(1)*b1, 0, Cn), 1, 2), 'LineColor', 'k', 'ShowText', 'Off');

%plot([d(1)*b1(1)], [d(1)*b1(2)], 'o', 'Color', 'r')
%plot([d(2)*b2(1)], [d(2)*b2(2)], 'o', 'Color', 'b')

set(gca, 'XTick', []);
set(gca, 'YTick', []);
set(gca, 'XLabel', []);
set(gca, 'YLabel', []);

print(handle, '-dpdf', '-r1500', fullfile(dp_path, 'cma_est_a_02'));


clf;
% objective contours
contour(X, Y, reshape(sqrt((X-opt(1)).^2+(Y-opt(2)).^2), length(X), length(Y)), 'LineColor', [0.75 0.75 0.75], 'LevelStep', 0.75);
hold on;
axis equal;

% update mean
mn = [0 0];
for i = 1:mu
  mn = mn + w(i)*pop(i,:);
end
plot([mn(1)], [mn(2)], '+', 'Color', 'black', 'MarkerSize', 6, 'MarkerFaceColor', 'k');

% population
plot(pop(:,1), pop(:,2), 'o', 'Color', 'black', 'MarkerFaceColor', 'k', 'MarkerSize', 1.5)

% CM contour
contour(X, Y, reshape(mvnpdf([X(:) Y(:)], 0, C), length(X), length(Y)), repmat(mvnpdf(sqrt([0 sigma]), 0, C), 1, 2), 'LineStyle', '--', 'LineColor', 'k', 'ShowText', 'Off');
contour(X, Y, reshape(mvnpdf([X(:) Y(:)], mn, Cn), length(X), length(Y)), repmat(mvnpdf(mn'+d(1)*b1, mn', Cn), 1, 2), 'LineColor', 'k', 'ShowText', 'Off');

set(gca, 'XTick', []);
set(gca, 'YTick', []);
set(gca, 'XLabel', []);
set(gca, 'YLabel', []);

print(handle, '-dpdf', '-r1500', fullfile(dp_path, 'cma_est_a_03'));