% S4 Fig.  Reported over calculated cascade coefficient against the ladder factor nu_2 of the
% twelve-site downstream ladder, at PII 36 (black), 5 (red) and 0.5 uM (blue); dashed line at 1.
S = plosstyle();
T = S.read('FigS4_ratio_vs_nu2.csv');
fig = figure('Units', 'centimeters', 'Position', [2 2 8.5 6.2], 'Color', 'w');
ax = axes('Position', [0.16 0.19 0.80 0.76]); hold(ax, 'on');
plot(ax, [1 max(T.nu2)], [1 1], '--', 'Color', S.gr, 'LineWidth', 1.0);
plot(ax, T.nu2, T.ratio_PII_36, '-', 'Color', S.k, 'LineWidth', 1.5, 'DisplayName', '36 \muM');
plot(ax, T.nu2, T.ratio_PII_5, '-', 'Color', S.r, 'LineWidth', 1.5, 'DisplayName', '5 \muM');
plot(ax, T.nu2, T.ratio_PII_0_5, '-', 'Color', S.b, 'LineWidth', 1.5, 'DisplayName', '0.5 \muM');
xlim(ax, [1 max(T.nu2)]); ylim(ax, [0 1.8]); yticks(ax, [0 0.5 1 1.5]);
xlabel(ax, '\nu_2'); ylabel(ax, 'reported / calculated');
legend(ax, 'Location', 'northeast', 'Box', 'on', 'EdgeColor', 'k'); S.axes(ax);
S.save(fig, 'FigS4');
