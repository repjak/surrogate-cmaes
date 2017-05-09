close all;

xl = 1;
density =1e5;
X = linspace(0, xl, density);
steepness = 2;
tr = 1;
fontSize = 20;
dp_path = '/home/jakub/Documents/dp/img';

handle = figure();
set(handle, 'PaperPositionMode', 'auto');

Xtrunc = min(X, tr)/tr;

y2 = arrayfun(@(x) GenerationsUpdater.simplesig(x, steepness), Xtrunc);
%y3 = arrayfun(@(x) GenerationsUpdater.simplesig((exp(x)-1)/(exp(1)-1), 2), 1-Xtrunc);
%y4 = arrayfun(@(x) (10^x-1), Xtrunc);

plot(X, X, '-k');
hold on;
plot(X, y2, '--k');
%plot(X, y3, '-.k');
%plot(X, y4, 'r');

axis equal;
xlim([0 xl]);
ylim([0 xl]);

set(gca, 'FontSize', fontSize);
set(gca, 'TickLabelInterpreter', 'latex');

xticks([0:0.2:1]);
yticks([0:0.2:1]);
hold on;
grid;
% xticklabels({'$0$', '$\frac{\varepsilon_T}{2}$', '$\varepsilon_T$'});
% yticklabels({'$0$', '$\frac{g_m^\mathrm{max}}{2}$', '$g_m^\mathrm{max}$'});
% xlabel('$\varepsilon$', 'Interpreter', 'latex');
% ylabel('$g_m$', 'Interpreter', 'latex');

lgd = legend('$T_1$', '$T_2,\,k=2$');
set(lgd, 'Interpreter', 'latex');

print(handle, '-depsc', fullfile(dp_path, 'transfer'));