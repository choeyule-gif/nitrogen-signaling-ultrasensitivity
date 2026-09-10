function Fig4()
%FIG4  Fig 4: the composite cascade response, one panel per file.
%   Fig4A  calibrated composite responses at PII 0.5, 5 and 36 uM with the
%          reported midpoints (circles)
%   Fig4B  calculated and reported coefficients at the three PII levels
%   Fig4C  reported over calculated, with the central 95% ranges of the
%          upstream coordinate perturbations
%
%   Reads figdata/Fig4*.csv and writes figures/Fig4A, Fig4B, Fig4C (PDF + PNG).
%
%   See also RUN_ALL_FIGURES, NEWPANEL, SAVEPANEL, SERIESCOLORS.

piiLevels = [0.5, 5, 36];    % uM
piiLabels = {'0.5', '5', '36'};
grey = [0.5, 0.5, 0.5];
thinEdge = 0.75;             % marker outlines, as in the reference figure

%% Fig4A: composite responses with the reported midpoints
composites = readFigData('Fig4A_composites.csv');
midpoints = readFigData('Fig4A_midpoints.csv');
colors = seriesColors(numel(piiLevels));

[fig, ax] = newPanel();
for k = 1:numel(piiLevels)
    plot(ax, composites.glutamine_mM, composites.(sprintf('Y_norm_PII_%g', piiLevels(k))), ...
        '-', 'Color', colors(k, :), 'DisplayName', [piiLabels{k}, ' \muM']);
end
for k = 1:numel(piiLevels)
    midpoint = midpoints.midpoint_mM(midpoints.PII_uM == piiLevels(k));
    plot(ax, midpoint, 0.5, 'o', 'Color', colors(k, :), 'MarkerFaceColor', 'w', ...
        'LineWidth', thinEdge, 'HandleVisibility', 'off');
end
ax.XScale = 'log';
xlim(ax, [1e-2, 10]);
ylim(ax, [-0.02, 1.02]);
xlabel(ax, '[Gln] / mM');
ylabel(ax, 'Y');
legend(ax, 'Location', 'northwest');
savePanel(fig, 'Fig4A');

%% Fig4B: calculated and reported coefficients
coefficients = readFigData('Fig4B_coefficients.csv');
colors = barColors();

[fig, ax] = newPanel();
bars = bar(ax, [coefficients.calculated, coefficients.reported], 'EdgeColor', 'none');
bars(1).FaceColor = colors(1, :);
bars(1).DisplayName = 'calculated';
bars(2).FaceColor = colors(2, :);
bars(2).DisplayName = 'reported';
xlim(ax, [0.4, 3.6]);
ylim(ax, [0, 8]);
xticks(ax, 1:numel(piiLevels));
xticklabels(ax, piiLabels);
xlabel(ax, 'T_{PII} / \muM');
ylabel(ax, 'n_H');
legend(ax, 'Location', 'northwest');
savePanel(fig, 'Fig4B');

%% Fig4C: reported over calculated, with the central 95% ranges
ratios = readFigData('Fig4C_ratios.csv');

[fig, ax] = newPanel();
errorbar(ax, ratios.PII_uM, ratios.ratio, ratios.ratio - ratios.p025, ...
    ratios.p975 - ratios.ratio, '-o', 'Color', seriesColors(1), 'MarkerFaceColor', 'w', ...
    'LineWidth', thinEdge, 'MarkerSize', 6, 'CapSize', 5);
yline(ax, 1, '--', 'Color', grey, 'LineWidth', 1);
ax.XScale = 'log';
xlim(ax, [0.2, 100]);
ylim(ax, [0.6, 2]);
xlabel(ax, 'T_{PII} / \muM');
ylabel(ax, 'reported / calculated');
savePanel(fig, 'Fig4C');
end
