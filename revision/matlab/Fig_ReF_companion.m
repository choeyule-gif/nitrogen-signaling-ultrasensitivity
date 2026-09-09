% Companion manuscript (stability): Re F(i omega) across n.  Blue, the range covered by the
% proofs (all omega at n = 1; omega above the circulation threshold otherwise); red, the band
% that remains open (n >= 2, omega below the threshold); circles at omega -> 0 mark F(0).
S = plosstyle();
T = S.read('Fig_ReF_companion.csv'); Mk = S.read('Fig_ReF_markers.csv');
ns = unique(T.n)';
fig = figure('Units', 'centimeters', 'Position', [2 2 11.5 8.5], 'Color', 'w');
ax = axes('Position', [0.02 0.05 0.96 0.92]); hold(ax, 'on');
for n = ns
    m = T.n == n; om = T.omega(m); z = T.ReF(m); lo = log10(om);
    thr = Mk.circulation_threshold(Mk.n == n); k = find(om >= thr, 1);
    plot3(ax, lo(k:end), n*ones(1, numel(lo)-k+1), z(k:end), '-', 'Color', S.b, 'LineWidth', 1.6);
    if n == 1
        plot3(ax, lo(1:k), n*ones(1, k), z(1:k), '-', 'Color', S.b, 'LineWidth', 1.6);
    else
        plot3(ax, lo(1:k), n*ones(1, k), z(1:k), '-', 'Color', S.r, 'LineWidth', 2.4);
    end
    plot3(ax, lo(1), n, Mk.F0(Mk.n == n), 'o', 'Color', S.b, 'MarkerSize', 4, 'MarkerFaceColor', 'w', 'LineWidth', 1.0);
end
plot3(ax, [-3 3], [0.3 0.3], [0 0], '-', 'Color', S.gr, 'LineWidth', 0.8);
xlabel(ax, 'log_{10}\omega'); ylabel(ax, 'n'); zlabel(ax, 'Re F');
xticks(ax, [-3 -1 1 3]); yticks(ax, ns); zticks(ax, [0 0.5 1]); zlim(ax, [0 1.05]);
set(ax, 'FontSize', 8, 'FontName', 'Helvetica', 'LineWidth', 0.8, 'Box', 'off'); grid(ax, 'off');
view(ax, -50, 20);
S.save(fig, 'Fig_ReF_companion');
