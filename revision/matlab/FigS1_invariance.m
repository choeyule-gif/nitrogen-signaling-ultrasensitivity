% S1 Fig.  (A) ratios s_i/s_{i-1} at n = 4 (non-identical ladder) against the free enzyme over
% a 5,000-fold range, solved in the raw species variables without imposing eq (1);
% (B) the same ratios at n = 3 while a dead-end complex E.S_0 has its association constant
% swept over eleven decades (T_E/T_L = 1.1e-9).
S = plosstyle();
fig = figure('Units', 'centimeters', 'Position', [2 2 17.4 6.2], 'Color', 'w');
cols = {S.k, S.r, S.b, S.g};
A = S.read('FigS1A_enzyme_sweep.csv');
ax = axes('Position', [0.08 0.19 0.38 0.75]); hold(ax, 'on');
for i = 1:4
    plot(ax, A.free_enzyme, A.(sprintf('s%d_over_s%d', i, i-1)), '-', 'Color', cols{i}, 'LineWidth', 1.5, 'DisplayName', sprintf('s_%d/s_%d', i, i-1));
end
set(ax, 'XScale', 'log', 'YScale', 'log'); xlim(ax, [min(A.free_enzyme) max(A.free_enzyme)]); ylim(ax, [0.1 10]);
yticks(ax, [0.1 1 10]); xlabel(ax, 'e'); ylabel(ax, 's_i/s_{i-1}');
legend(ax, 'Location', 'south', 'NumColumns', 2, 'Box', 'on', 'EdgeColor', 'k'); S.axes(ax); S.label(ax, 'A');
B = S.read('FigS1B_deadend_sweep.csv');
ax = axes('Position', [0.59 0.19 0.38 0.75]); hold(ax, 'on');
for i = 1:3
    plot(ax, B.K_Z, B.(sprintf('s%d_over_s%d', i, i-1)), '-', 'Color', cols{i}, 'LineWidth', 1.5, 'DisplayName', sprintf('s_%d/s_%d', i, i-1));
end
set(ax, 'XScale', 'log', 'YScale', 'log'); xlim(ax, [1e-6 1e5]); ylim(ax, [0.25 5]);
xticks(ax, [1e-6 1e-3 1 1e3]); yticks(ax, [0.3 1 3]); xlabel(ax, 'K_Z'); ylabel(ax, 's_i/s_{i-1}');
legend(ax, 'Location', 'east', 'Box', 'on', 'EdgeColor', 'k'); S.axes(ax); S.label(ax, 'B');
S.save(fig, 'FigS1');
