% Fig 1.  (A) the two cycles of the cascade (schematic); (B) mean modified fraction of the
% free target for n = 1, 2, 3, 12 identical independent sites, which coincide; (C) the local
% coefficient of the same curves, equal to 1.
S = plosstyle();
fig = figure('Units', 'centimeters', 'Position', [2 2 17.4 9.2], 'Color', 'w');

% ---- (A) schematic: boxes and arrows, labels only
ax = axes('Position', [0.05 0.62 0.90 0.33]); axis(ax, [0 1 0 1]); axis(ax, 'off'); hold(ax, 'on');
S.label(ax, 'A'); set(findobj(ax, 'Type', 'text'), 'Position', [0.0 1.02 0]);
xs = [0.02 0.235 0.45 0.665 0.88]; w = 0.12; h = 0.42; y0 = 0.30;
lab = {{'Gln'}, {'GlnD', 'UT / UR'}, {'PII', '0–3 UMP'}, {'ATase', 'AT / AR'}, {'GS', '0–12 AMP'}};
face = {[1 1 1], S.lb, [1 1 1], S.lb, [1 1 1]};
for i = 1:5
    rectangle(ax, 'Position', [xs(i) y0 w h], 'Curvature', [0.25 0.5], 'FaceColor', face{i}, 'EdgeColor', [0.3 0.3 0.3], 'LineWidth', 0.9);
    text(ax, xs(i)+w/2, y0+h/2, lab{i}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 9, 'FontName', 'Helvetica');
end
for i = 1:4        % arrows drawn as line + triangle in data units
    xa = xs(i)+w+0.005; xb = xs(i+1)-0.005; ym = y0+h/2;
    plot(ax, [xa xb-0.02], [ym ym], '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.3);
    patch(ax, [xb-0.025 xb xb-0.025], [ym-0.06 ym ym+0.06], [0.3 0.3 0.3], 'EdgeColor', 'none');
end
% under-labels: the two factors and the collective readout (symbols only)
text(ax, xs(2)+w/2, 0.10, '\eta', 'HorizontalAlignment', 'center', 'FontSize', 11, 'Color', S.b, 'FontName', 'Helvetica');
text(ax, xs(3)+w/2, 0.10, '\nu', 'HorizontalAlignment', 'center', 'FontSize', 11, 'Color', S.o, 'FontName', 'Helvetica');
text(ax, xs(4)+w/2, 0.10, 'p_0 = (1-\theta)^3', 'HorizontalAlignment', 'center', 'FontSize', 9, 'Color', S.g, 'FontName', 'Helvetica');

% ---- (B) theta(l) for n = 1, 2, 3, 12
T = S.read('Fig1B_theta.csv');
ax = axes('Position', [0.10 0.13 0.36 0.40]); hold(ax, 'on');
ns = [1 2 3 12]; cols = {S.k, S.r, S.b, S.g}; lw = [2.6 2.0 1.4 0.9];
for i = 1:4
    plot(ax, T.l, T.(sprintf('theta_n%d', ns(i))), '-', 'Color', cols{i}, 'LineWidth', lw(i), 'DisplayName', sprintf('n = %d', ns(i)));
end
set(ax, 'XScale', 'log'); xlim(ax, [10^-2.5 10^2.5]); ylim(ax, [-0.02 1.02]);
xticks(ax, [1e-2 1e0 1e2]); yticks(ax, [0 0.5 1]);
xlabel(ax, 'l'); ylabel(ax, '\theta'); legend(ax, 'Location', 'northeast', 'Box', 'on', 'EdgeColor', 'k');
S.axes(ax); S.label(ax, 'B');

% ---- (C) local coefficient
T = S.read('Fig1C_local_hill.csv');
ax = axes('Position', [0.60 0.13 0.36 0.40]); hold(ax, 'on');
for i = 1:4
    plot(ax, T.l, T.(sprintf('nH_n%d', ns(i))), '-', 'Color', cols{i}, 'LineWidth', lw(i));
end
set(ax, 'XScale', 'log'); xlim(ax, [10^-2.5 10^2.5]); ylim(ax, [0.98 1.02]);
xticks(ax, [1e-2 1e0 1e2]); yticks(ax, [0.98 1.00 1.02]);
xlabel(ax, 'l'); ylabel(ax, 'n_H^{loc}');
S.axes(ax); S.label(ax, 'C');
S.save(fig, 'Fig1');
