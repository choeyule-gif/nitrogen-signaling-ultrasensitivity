function Fig1()
%FIG1  Fig 1: n identical independent sites, one panel per file.
%   Fig1A  mean modified fraction theta of the free target against the input l
%          for n = 1, 2, 3 and 12 sites; the four curves coincide
%   Fig1B  local coefficient of the same curves, equal to 1 throughout
%
%   Reads figdata/Fig1B_theta.csv and figdata/Fig1C_local_hill.csv (the file
%   names keep the letters of the earlier layout, whose schematic panel A was
%   dropped) and writes figures/Fig1A and figures/Fig1B (PDF + PNG).
%
%   See also RUN_ALL_FIGURES, NEWPANEL, SAVEPANEL, SERIESCOLORS.

siteCounts = [1, 2, 3, 12];
lineStyles = {'-', '--', ':', '-.'};    % the curves coincide; styles tell them apart
colors = seriesColors(numel(siteCounts));
inputLimits = [10^-2.5, 10^2.5];

%% Fig1A: modified fraction against the input
theta = readFigData('Fig1B_theta.csv');

[fig, ax] = newPanel();
for k = 1:numel(siteCounts)
    plot(ax, theta.l, theta.(sprintf('theta_n%d', siteCounts(k))), lineStyles{k}, ...
        'Color', colors(k, :), 'DisplayName', sprintf('n = %d', siteCounts(k)));
end
ax.XScale = 'log';
xlim(ax, inputLimits);
ylim(ax, [-0.02, 1.02]);
xticks(ax, [1e-2, 1, 1e2]);
xlabel(ax, 'l');
ylabel(ax, '\theta');
legend(ax, 'Location', 'northeast');
savePanel(fig, 'Fig1A');

%% Fig1B: local coefficient of the same curves
localHill = readFigData('Fig1C_local_hill.csv');

[fig, ax] = newPanel();
for k = 1:numel(siteCounts)
    plot(ax, localHill.l, localHill.(sprintf('nH_n%d', siteCounts(k))), lineStyles{k}, ...
        'Color', colors(k, :), 'DisplayName', sprintf('n = %d', siteCounts(k)));
end
ax.XScale = 'log';
xlim(ax, inputLimits);
ylim(ax, [0.98, 1.02]);
xticks(ax, [1e-2, 1, 1e2]);
xlabel(ax, 'l');
ylabel(ax, 'n_H^{loc}');
legend(ax, 'Location', 'northeast');
savePanel(fig, 'Fig1B');
end
