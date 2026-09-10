function FigS4_nu2()
%FIGS4_NU2  S4 Fig: reported over calculated against the downstream ladder factor.
%   Reported over calculated cascade coefficient against the ladder factor nu_2
%   of the twelve-site downstream ladder at PII 36, 5 and 0.5 uM; dashed line
%   at 1.
%
%   Reads figdata/FigS4_ratio_vs_nu2.csv and writes figures/FigS4 (PDF + PNG).
%
%   See also RUN_ALL_FIGURES, NEWPANEL, SAVEPANEL, SERIESCOLORS.

ratioTable = readFigData('FigS4_ratio_vs_nu2.csv');
piiColumns = {'ratio_PII_36', 'ratio_PII_5', 'ratio_PII_0.5'};
piiLabels = {'36 \muM', '5 \muM', '0.5 \muM'};
colors = seriesColors(numel(piiColumns));
grey = [0.5, 0.5, 0.5];

[fig, ax] = newPanel();
for k = 1:numel(piiColumns)
    plot(ax, ratioTable.nu2, ratioTable.(piiColumns{k}), '-', 'Color', colors(k, :), ...
        'DisplayName', piiLabels{k});
end
yline(ax, 1, '--', 'Color', grey, 'LineWidth', 1, 'HandleVisibility', 'off');
xlim(ax, [1, max(ratioTable.nu2)]);
ylim(ax, [0, 1.8]);
xlabel(ax, '\nu_2');
ylabel(ax, 'reported / calculated');
legend(ax, 'Location', 'northeast');
savePanel(fig, 'FigS4');
end
