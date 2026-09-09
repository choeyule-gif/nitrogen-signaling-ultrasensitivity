% Fig 4.  (A) calibrated composite responses at PII 0.5, 5 and 36 uM with the reported
% midpoints; (B) calculated (filled) and reported (open) coefficients; (C) reported over
% calculated with the central 95% ranges of the upstream coordinate perturbations.
S = plosstyle();
fig = figure('Units', 'centimeters', 'Position', [2 2 17.4 5.9], 'Color', 'w');
pos = {[0.075 0.19 0.24 0.74], [0.415 0.19 0.24 0.74], [0.755 0.19 0.24 0.74]};
cols = {S.b, S.r, S.g}; pii = [0.5 5 36];

% ---- (A)
T = S.read('Fig4A_composites.csv'); Md = S.read('Fig4A_midpoints.csv');
ax = axes('Position', pos{1}); hold(ax, 'on');
for i = 1:3
    plot(ax, T.glutamine_mM, T.(sprintf('Y_norm_PII_%g', pii(i))), '-', 'Color', cols{i}, 'LineWidth', 1.5, 'DisplayName', sprintf('%g \\muM', pii(i)));
    plot(ax, Md.midpoint_mM(i), 0.5, 'o', 'Color', cols{i}, 'MarkerSize', 5, 'MarkerFaceColor', 'w', 'LineWidth', 1.1, 'HandleVisibility', 'off');
end
set(ax, 'XScale', 'log'); xlim(ax, [1e-2 10]); ylim(ax, [-0.02 1.02]);
xticks(ax, [1e-2 1e-1 1 10]); yticks(ax, [0 0.5 1]);
xlabel(ax, '[Gln] / mM'); ylabel(ax, 'Y'); legend(ax, 'Location', 'northwest', 'Box', 'on', 'EdgeColor', 'k');
S.axes(ax); S.label(ax, 'A');

% ---- (B)
C = S.read('Fig4B_coefficients.csv');
ax = axes('Position', pos{2}); hold(ax, 'on');
x = 1:3; w = 0.36;
bar(ax, x-w/2, C.calculated, w, 'FaceColor', S.b, 'EdgeColor', S.k, 'LineWidth', 0.7);
bar(ax, x+w/2, C.reported, w, 'FaceColor', 'w', 'EdgeColor', S.k, 'LineWidth', 1.0);
xlim(ax, [0.4 3.6]); ylim(ax, [0 8]); xticks(ax, x); xticklabels(ax, {'0.5', '5', '36'}); yticks(ax, [0 2 4 6 8]);
xlabel(ax, 'T_{PII} / \muM'); ylabel(ax, 'n_H'); S.axes(ax); set(ax, 'XMinorTick', 'off'); S.label(ax, 'B');

% ---- (C)
Rt = S.read('Fig4C_ratios.csv');
ax = axes('Position', pos{3}); hold(ax, 'on');
plot(ax, [0.2 100], [1 1], '--', 'Color', S.gr, 'LineWidth', 1.0);
errorbar(ax, Rt.PII_uM, Rt.ratio, Rt.ratio-Rt.p025, Rt.p975-Rt.ratio, 'o-', 'Color', S.r, 'LineWidth', 1.2, 'MarkerSize', 5, 'MarkerFaceColor', 'w', 'CapSize', 4);
set(ax, 'XScale', 'log'); xlim(ax, [0.2 100]); ylim(ax, [0.6 2.0]);
xticks(ax, [1 10 100]); yticks(ax, [0.8 1.2 1.6 2.0]);
xlabel(ax, 'T_{PII} / \muM'); ylabel(ax, 'reported / calculated'); S.axes(ax); S.label(ax, 'C');
S.save(fig, 'Fig4');
