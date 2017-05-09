%% init
phi = (1 + sqrt(5))/2;
ylim = 5;
[X, Y] = meshgrid(linspace(-(ylim*phi/2), phi*ylim/2, 1e2), linspace(-ylim/2, ylim/2, 1e2));
opt = [5*phi 5];
fitfun = @(x, y) sqrt((x-opt(1)).^2+(y-opt(2)).^2);
q = @(x) norm(x-opt)^2;
sigma = 1.2;
lambda = 12;
mu = ceil(0.25 * lambda);
crossSize = 15;
pointSize = 17;
textSize = 30;

%% step 1
angle = pi/3;
m1 = [-1.2, -0.5];
B = [cos(angle) sin(angle); -sin(angle) cos(angle)];
D = [0.75 0; 0 1];
C1 = B*D*D*B';

w = [1/(mu) 1/(4*mu) repmat((1-1/(mu)-1/(4*mu))/(mu-2), 1, mu-2)];

handle = figure('Units', 'inches');

dp_path = '/home/jakub/Documents/itat_slides/img/';

% objective contours
contour(X, Y, reshape(sqrt((X-opt(1)).^2+(Y-opt(2)).^2), length(X), length(Y)));
hold on;
axis equal;

% CM contour
contour(X, Y, reshape(mvnpdf([X(:) Y(:)], m1, C1), length(X), length(Y)), repmat(mvnpdf(m1'+D(1,1)*B(:,1), m1', C1), 1, 2), 'LineColor', 'k');
cax = m1'+D(1,1)*B(:,1);
h_C1_text = text(cax(1), cax(2), '$\sigma_1\mathbf{C}_1$', 'VerticalAlignment', 'top', 'FontSize', textSize, 'Interpreter', 'latex');

% mean
h_mean1 = plot(m1(1), m1(2), '+', 'Color', 'black', 'MarkerSize', crossSize, 'MarkerFaceColor', 'k');
h_mean1_text = text(m1(1), m1(2), '$\mathbf{m}_1$', 'VerticalAlignment', 'top', 'FontSize', textSize, 'Interpreter', 'latex');

% sampling
pop = mvnrnd(m1, sigma*C1, lambda);
plot(pop(:,1), pop(:,2), 'o', 'Color', 'b', 'MarkerSize', pointSize);

delete(h_mean1);
delete(h_mean1_text);

h_mean1 = plot(m1(1), m1(2), '+', 'Color', 'black', 'MarkerSize', crossSize, 'MarkerFaceColor', 'k');
h_mean1_text = text(m1(1), m1(2), '$\mathbf{m}_1$', 'VerticalAlignment', 'top', 'FontSize', textSize, 'Interpreter', 'latex');

set(gca, 'XTick', []);
set(gca, 'YTick', []);
set(gca, 'XLabel', []);
set(gca, 'YLabel', []);

set(handle,'PaperPositionMode','manual');
set(handle, 'PaperUnits','inches');
set(handle, 'Units','inches');
pos=get(handle,'Position');
set(handle, 'PaperSize', [5*phi 5]);
set(handle, 'PaperPosition', [0 0 5*phi 5]);

print(handle, '-dpdf', '-r1500', fullfile(dp_path, 'cma_01'));


%% evaluation
cmap = colormap;
for i=1:lambda
  color_ind = round(size(colormap, 1) * min(1, fitfun(pop(i, 1), pop(i, 2)) / fitfun(-ylim*phi/2, -ylim/2)));
  pop_colors_rgb(i, :) = cmap(color_ind, :);
end

for i=1:lambda
  plot(pop(i,1), pop(i,2), 'o', 'Color', 'b', 'MarkerSize', pointSize, 'MarkerFaceColor', pop_colors_rgb(i, :));
end

delete(h_mean1);
delete(h_mean1_text);
delete(h_C1_text);

h_mean1 = plot(m1(1), m1(2), '+', 'Color', 'black', 'MarkerSize', crossSize, 'MarkerFaceColor', 'k');
h_mean1_text = text(m1(1), m1(2), '$\mathbf{m}_1$', 'VerticalAlignment', 'top', 'FontSize', textSize, 'Interpreter', 'latex');

contour(X, Y, reshape(mvnpdf([X(:) Y(:)], m1, C1), length(X), length(Y)), repmat(mvnpdf(m1'+D(1,1)*B(:,1), m1', C1), 1, 2), 'LineColor', 'k');
cax = m1'+D(1,1)*B(:,1);
h_C1_text = text(cax(1), cax(2), '$\sigma_1\mathbf{C}_1$', 'VerticalAlignment', 'top', 'FontSize', textSize, 'Interpreter', 'latex');

print(handle, '-dpdf', '-r1500', fullfile(dp_path, 'cma_02'));


%% selection
f = [];
for x = pop'
  f(end+1) = q(x');
end

[s, i] = sort(f);
pop = pop(i(1:mu),:);
pop_colors_rgb = pop_colors_rgb(i(1:mu),:);

for i=1:mu
  plot(pop(i,1), pop(i,2), 'o', 'Color', 'r', 'MarkerSize', pointSize, 'MarkerFaceColor', pop_colors_rgb(i, :));
end

delete(h_mean1);
delete(h_mean1_text);

h_mean1 = plot(m1(1), m1(2), '+', 'Color', 'black', 'MarkerSize', crossSize, 'MarkerFaceColor', 'k');
h_mean1_text = text(m1(1), m1(2), '$\mathbf{m}_1$', 'VerticalAlignment', 'top', 'FontSize', textSize, 'Interpreter', 'latex');

contour(X, Y, reshape(mvnpdf([X(:) Y(:)], m1, C1), length(X), length(Y)), repmat(mvnpdf(m1'+D(1,1)*B(:,1), m1', C1), 1, 2), 'LineColor', 'k');
cax = m1'+D(1,1)*B(:,1);
h_C1_text = text(cax(1), cax(2), '$\sigma_1\mathbf{C}_1$', 'VerticalAlignment', 'top', 'FontSize', textSize, 'Interpreter', 'latex');

print(handle, '-dpdf', '-r1500', fullfile(dp_path, 'cma_03'));

%% update mean
% for x = pop'
%   y = [m1(1) x(1)];
%   z = [m1(2) x(2)];
%   %quiver(z(1), y(1), z(2)-z(1), y(2)-y(1), 0, 'Color', 'k');
%   arrow(y, z);
%   hold on;
% end

delete(h_mean1);
delete(h_mean1_text);

plot([m1(1)], [m1(2)], '+', 'Color', 0.3*ones(1,3), 'MarkerSize', crossSize, 'MarkerFaceColor', 'k');
text(m1(1), m1(2), '$\mathbf{m}_1$', 'VerticalAlignment', 'top', 'FontSize', textSize, 'Interpreter', 'latex', 'Color', 0.3*ones(1,3));

m2 = [0 0];
for i = 1:mu
  m2 = m2 + w(i)*pop(i,:);
end
plot([m2(1)], [m2(2)], '+', 'Color', 'black', 'MarkerSize', crossSize, 'MarkerFaceColor', 'k');
text(m2(1), m2(2), '$\mathbf{m}_2$', 'VerticalAlignment', 'top', 'FontSize', textSize, 'Interpreter', 'latex');

print(handle, '-dpdf', '-r1500', fullfile(dp_path, 'cma_04'));

%% update C
%delete(h_C1);
delete(h_C1_text);

contour(X, Y, reshape(mvnpdf([X(:) Y(:)], m1, C1), length(X), length(Y)), repmat(mvnpdf(m1'+D(1,1)*B(:,1), m1', C1), 1, 2), 'LineColor', 0.3*ones(1,3));
cax = m1'+D(1,1)*B(:,1);
text(cax(1), cax(2), '$\sigma_1\mathbf{C}_1$', 'VerticalAlignment', 'top', 'FontSize', textSize, 'Interpreter', 'latex', 'Color', 0.3*ones(1,3));

C2 = zeros(2, 2);
for i = 1:mu
  C2 = C2 + w(i)*(pop(i,:)-m1)'*(pop(i,:)-m1);
end
%Cn = (1/(sigma^2))*Cn;

[B, D] = eig(C2);

% CM contour
contour(X, Y, reshape(mvnpdf([X(:) Y(:)], m2, C2), length(X), length(Y)), repmat(mvnpdf(m2'+sqrt(D(1,1))*B(:,1), m2', C2), 1, 2), 'LineColor', 'k', 'ShowText', 'Off');
cax = m2'-sqrt(D(1,1))*B(:,1);
text(cax(1), cax(2), '$\sigma_2\mathbf{C}_2$', 'VerticalAlignment', 'bottom', 'FontSize', textSize, 'Interpreter', 'latex');

print(handle, '-dpdf', '-r1500', fullfile(dp_path, 'cma_05'));

