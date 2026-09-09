% Fig 2.  (A) curves of constant local coefficient nu*eta in the (nu, eta) plane with the
% independent ladder nu = 1 (dashed) and the trimer ceiling nu = 3; (B) the maximum switch
% elasticity of eq (6) against the affinity ratio, with the Mg (circle) and Mn (square)
% ratios; (C) the local logit slope of the fitted first-layer curve against glutamine, the
% fitted exponent (dashed) and the value at the midpoint (circle).
S = plosstyle();
fig = figure('Units', 'centimeters', 'Position', [2 2 17.4 5.9], 'Color', 'w');
pos = {[0.075 0.19 0.24 0.74], [0.415 0.19 0.24 0.74], [0.755 0.19 0.24 0.74]};

% ---- (A)
T = S.read('Fig2A_isocurves.csv');
ax = axes('Position', pos{1}); hold(ax, 'on');
cs = [1 2 3 4 6];
for c = cs
    e = T.(sprintf('eta_c%d', c)); m = e >= 0.5 & e <= 3.5;
    plot(ax, T.nu(m), e(m), '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.8);
    xm = min(2.8, c/0.65); ym = c/xm;
    text(ax, xm, ym, sprintf('%d', c), 'FontSize', 7, 'Color', S.gr, 'HorizontalAlignment', 'center', 'BackgroundColor', 'w', 'Margin', 0.5);
end
plot(ax, [1 1], [0.5 3.5], '--', 'Color', S.b, 'LineWidth', 1.0);
plot(ax, [3 3], [0.5 3.5], '-', 'Color', S.p, 'LineWidth', 1.4);
xlim(ax, [0.5 3.1]); ylim(ax, [0.5 3.5]); xticks(ax, [1 2 3]); yticks(ax, [1 2 3]);
xlabel(ax, '\nu'); ylabel(ax, '\eta'); S.axes(ax); S.label(ax, 'A');

% ---- (B)
T = S.read('Fig2B_eta_max.csv'); Mk = S.read('Fig2B_markers.csv');
ax = axes('Position', pos{2}); hold(ax, 'on');
plot(ax, [1 1e4], [2 2], ':', 'Color', S.gr, 'LineWidth', 0.9);
plot(ax, T.affinity_ratio, T.eta_max, '-', 'Color', S.r, 'LineWidth', 1.5);
plot(ax, Mk.affinity_ratio(1), Mk.eta_max(1), 'o', 'Color', S.k, 'MarkerSize', 5, 'MarkerFaceColor', 'w', 'LineWidth', 1.1);
plot(ax, Mk.affinity_ratio(2), Mk.eta_max(2), 's', 'Color', S.k, 'MarkerSize', 5.5, 'MarkerFaceColor', 'w', 'LineWidth', 1.1);
set(ax, 'XScale', 'log'); xlim(ax, [1 1e4]); ylim(ax, [0.98 2.04]);
xticks(ax, [1 1e2 1e4]); yticks(ax, [1 1.5 2]);
xlabel(ax, 'A_{UT}/A_{UR}'); ylabel(ax, 'max \eta'); S.axes(ax); S.label(ax, 'B');

% ---- (C)
T = S.read('Fig2C_local_slope.csv'); Mk = S.read('Fig2C_markers.csv');
Sm = Mk.value(1); nloc = Mk.value(2); hfit = Mk.value(3);
ax = axes('Position', pos{3}); hold(ax, 'on');
plot(ax, [1e-2 10^1.3], [hfit hfit], '--', 'Color', S.r, 'LineWidth', 1.0);
plot(ax, T.glutamine_mM, T.nH_loc, '-', 'Color', S.b, 'LineWidth', 1.5);
plot(ax, Sm, nloc, 'o', 'Color', S.b, 'MarkerSize', 5, 'MarkerFaceColor', 'w', 'LineWidth', 1.1);
set(ax, 'XScale', 'log'); xlim(ax, [1e-2 10^1.3]); ylim(ax, [0 2.3]);
xticks(ax, [1e-2 1e-1 1 10]); yticks(ax, [0 1 2]);
xlabel(ax, '[Gln] / mM'); ylabel(ax, 'n_H'); S.axes(ax); S.label(ax, 'C');
S.save(fig, 'Fig2');
