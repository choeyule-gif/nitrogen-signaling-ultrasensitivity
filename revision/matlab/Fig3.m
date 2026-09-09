% Fig 3.  (A) the digitised titration, the bounded four-parameter fit and the pointwise
% 95% band under the coordinate perturbations; (B) standardised residuals of the fixed-
% exponent fits; (C) the fitted exponent over the 2,000 perturbations; (D) profile
% likelihood loss against the regulatory swing, with the criterion and the lower endpoint.
S = plosstyle();
fig = figure('Units', 'centimeters', 'Position', [2 2 17.4 12.0], 'Color', 'w');
pos = {[0.085 0.60 0.38 0.35], [0.585 0.60 0.38 0.35], [0.085 0.10 0.38 0.35], [0.585 0.10 0.38 0.35]};

% ---- (A)
P = S.read('Fig3A_points.csv'); F = S.read('Fig3A_fit.csv');
ax = axes('Position', pos{1}); hold(ax, 'on');
fill(ax, [F.glutamine_mM; flipud(F.glutamine_mM)], [F.band_lo; flipud(F.band_hi)], S.lb, 'EdgeColor', 'none');
plot(ax, F.glutamine_mM, F.fit, '-', 'Color', S.b, 'LineWidth', 1.5);
plot(ax, P.glutamine_mM, P.UMP_per_trimer, 'o', 'Color', S.k, 'MarkerSize', 5, 'MarkerFaceColor', 'w', 'LineWidth', 1.1);
set(ax, 'XScale', 'log'); xlim(ax, [1e-2 10^1.25]); ylim(ax, [0 3.15]);
xticks(ax, [1e-2 1e-1 1 10]); yticks(ax, [0 1 2 3]);
xlabel(ax, '[Gln] / mM'); ylabel(ax, 'UMP per trimer'); S.axes(ax); S.label(ax, 'A');

% ---- (B)
R = S.read('Fig3B_residuals.csv');
ax = axes('Position', pos{2}); hold(ax, 'on');
fill(ax, [1e-2 20 20 1e-2], [-2 -2 2 2], [0.9 0.9 0.9], 'EdgeColor', 'none');
plot(ax, [1e-2 20], [0 0], '-', 'Color', S.k, 'LineWidth', 0.7);
hs = [1 2 3]; cols = {S.gr, S.b, S.r};
for i = 1:3
    plot(ax, R.glutamine_mM, R.(sprintf('r_h%d', hs(i))), 'o-', 'Color', cols{i}, 'LineWidth', 1.0, 'MarkerSize', 4, 'MarkerFaceColor', 'w', 'DisplayName', sprintf('h = %d', hs(i)));
end
set(ax, 'XScale', 'log'); xlim(ax, [1e-2 20]); ylim(ax, [-8 11]);
xticks(ax, [1e-2 1e-1 1 10]); yticks(ax, [-8 -4 0 4 8]);
xlabel(ax, '[Gln] / mM'); ylabel(ax, 'r / \sigma'); legend(ax, 'Location', 'northeast', 'Box', 'on', 'EdgeColor', 'k', 'NumColumns', 3);
S.axes(ax); S.label(ax, 'B');

% ---- (C)
H = S.read('Fig3C_exponent_draws.csv'); Mk = S.read('Fig3C_markers.csv');
ax = axes('Position', pos{3}); hold(ax, 'on');
histogram(ax, H{:,1}, 35, 'Normalization', 'pdf', 'FaceColor', S.lb, 'EdgeColor', 'w');
plot(ax, [Mk.value(1) Mk.value(1)], [0 3], '-', 'Color', S.b, 'LineWidth', 1.3);
plot(ax, [Mk.value(2) Mk.value(2)], [0 3], ':', 'Color', S.gr, 'LineWidth', 1.1);
plot(ax, [Mk.value(3) Mk.value(3)], [0 3], ':', 'Color', S.gr, 'LineWidth', 1.1);
xlim(ax, [1.5 3.0]); ylim(ax, [0 2.8]); xticks(ax, [1.5 2 2.5 3]); yticks(ax, [0 1 2]);
xlabel(ax, 'h'); ylabel(ax, 'density'); S.axes(ax); S.label(ax, 'C');

% ---- (D)
W = S.read('Fig3D_swing_profile.csv'); Mk = S.read('Fig3D_markers.csv');
ax = axes('Position', pos{4}); hold(ax, 'on');
plot(ax, [10^1.5 1e5], [Mk.value(1) Mk.value(1)], '--', 'Color', S.gr, 'LineWidth', 1.0);
plot(ax, [Mk.value(2) Mk.value(2)], [0 12], ':', 'Color', S.p, 'LineWidth', 1.1);
plot(ax, W.swing, W.loss, '-', 'Color', S.p, 'LineWidth', 1.5);
set(ax, 'XScale', 'log'); xlim(ax, [10^1.5 1e5]); ylim(ax, [0 12]);
xticks(ax, [1e2 1e3 1e4 1e5]); yticks(ax, [0 4 8 12]);
xlabel(ax, 'R_{fwd} R_{rev}'); ylabel(ax, '2\Delta\itL'); S.axes(ax); S.label(ax, 'D');
S.save(fig, 'Fig3');
