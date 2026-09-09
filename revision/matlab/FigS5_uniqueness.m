% S5 Fig.  The residual R(l) whose zeros are the positive steady states: (A) both conditions
% of the uniqueness criterion satisfied (one zero, R increasing); (B) condition (i) broken by
% a spread of 2.89 in the catalytic ratios at n = 4 (three zeros).
S = plosstyle();
fig = figure('Units', 'centimeters', 'Position', [2 2 17.4 6.2], 'Color', 'w');
A = S.read('FigS5A_residual.csv');
ax = axes('Position', [0.08 0.19 0.38 0.75]); hold(ax, 'on');
plot(ax, [1e-2 1e2], [0 0], '--', 'Color', S.k, 'LineWidth', 0.8);
plot(ax, A.l, A.R, '-', 'Color', S.b, 'LineWidth', 1.5);
k = find(diff(sign(A.R)) ~= 0); plot(ax, A.l(k), zeros(size(k)), 'o', 'Color', S.r, 'MarkerSize', 6, 'MarkerFaceColor', S.r);
set(ax, 'XScale', 'log'); xlim(ax, [1e-2 1e2]); ylim(ax, [-1.5 6]); xticks(ax, [1e-2 1 1e2]);
xlabel(ax, 'l'); ylabel(ax, 'R(l)'); S.axes(ax); S.label(ax, 'A');
B = S.read('FigS5B_residual.csv');
ax = axes('Position', [0.59 0.19 0.38 0.75]); hold(ax, 'on');
plot(ax, [min(B.l) max(B.l)], [0 0], '--', 'Color', S.k, 'LineWidth', 0.8);
plot(ax, B.l, B.R, '-', 'Color', S.b, 'LineWidth', 1.5);
k = find(diff(sign(B.R)) ~= 0); plot(ax, B.l(k), zeros(size(k)), 'o', 'Color', S.r, 'MarkerSize', 6, 'MarkerFaceColor', S.r);
set(ax, 'XScale', 'log'); xlim(ax, [min(B.l) max(B.l)]); xticks(ax, [1e-4 1e-2 1]);
xlabel(ax, 'l'); ylabel(ax, 'R(l)'); S.axes(ax); S.label(ax, 'B');
S.save(fig, 'FigS5');
