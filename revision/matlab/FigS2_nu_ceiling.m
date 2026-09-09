% S2 Fig.  Supremum over the input of the ladder factor nu for every three-site ladder
% c = (1, 3a, 3b, 1), as a surface over (log10 a, log10 b); a = b = 1 is the binomial ladder
% (nu = 1, circle); the ceiling nu = 3 is approached as a, b -> 0.  The two lines are the
% ladders reproducing the observed ratio 1.1106 with the fitted-exponent (blue) and range
% (red) estimators.  Rotate interactively and re-export.
S = plosstyle();
T = S.read('FigS2_nu_sup_surface.csv');
la = unique(T.log10_a); lb = unique(T.log10_b); na = numel(la); nb = numel(lb);
Z = reshape(T.nu_sup, na, nb)';                    % rows: log10 b, columns: log10 a
fig = figure('Units', 'centimeters', 'Position', [2 2 11.5 8.5], 'Color', 'w');
ax = axes('Position', [0.02 0.05 0.96 0.92]); hold(ax, 'on');
surf(ax, la, lb, Z, 'EdgeColor', [0.35 0.35 0.35], 'LineWidth', 0.25);
colormap(ax, parula); clim(ax, [1 3]);
plot3(ax, 0, 0, 1, 'o', 'Color', S.k, 'MarkerSize', 6, 'MarkerFaceColor', 'w', 'LineWidth', 1.2);
Lf = S.read('FigS2_locus_fit.csv'); Lr = S.read('FigS2_locus_range.csv');
plot3(ax, Lf.log10_a, Lf.log10_b, Lf.nu_sup+0.02, '-', 'Color', S.b, 'LineWidth', 1.8);
plot3(ax, Lr.log10_a, Lr.log10_b, Lr.nu_sup+0.02, '-', 'Color', S.r, 'LineWidth', 1.8);
xlabel(ax, 'log_{10} a'); ylabel(ax, 'log_{10} b'); zlabel(ax, 'sup_u \nu');
xticks(ax, [-3 0 3]); yticks(ax, [-3 0 3]); zticks(ax, [1 2 3]); zlim(ax, [1 3.05]);
set(ax, 'FontSize', 8, 'FontName', 'Helvetica', 'LineWidth', 0.8, 'Box', 'off'); grid(ax, 'off');
view(ax, -40, 24);
S.save(fig, 'FigS2');
