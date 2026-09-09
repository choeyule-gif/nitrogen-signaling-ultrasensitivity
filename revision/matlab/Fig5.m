% Fig 5.  (A) ligands bound per enzyme for a one-site and a concerted two-site binding law;
% (B) glutamine titrations predicted by the minimal network at 1x/2x enzyme and 1x/4x
% target (the four curves coincide); (C) the same conditions with a catalytically active
% complex bridging two targets; (D) twelve-site distributions with the same mean for
% nu = 1 and nu = 2.2.
S = plosstyle();
fig = figure('Units', 'centimeters', 'Position', [2 2 17.4 11.4], 'Color', 'w');
pos = {[0.085 0.60 0.38 0.35], [0.585 0.60 0.38 0.35], [0.085 0.10 0.38 0.35], [0.585 0.10 0.38 0.35]};

% ---- (A)
T = S.read('Fig5A_binding.csv');
ax = axes('Position', pos{1}); hold(ax, 'on');
plot(ax, T.scaled_glutamine, T.one_ligand, '-', 'Color', S.b, 'LineWidth', 1.5, 'DisplayName', 'h = 1');
plot(ax, T.scaled_glutamine, T.two_ligands, '-', 'Color', S.r, 'LineWidth', 1.5, 'DisplayName', 'h = 2');
set(ax, 'XScale', 'log'); xlim(ax, [1e-2 1e2]); ylim(ax, [0 2.05]);
xticks(ax, [1e-2 1 1e2]); yticks(ax, [0 1 2]);
xlabel(ax, '[Gln] / K_L'); ylabel(ax, 'bound / enzyme'); legend(ax, 'Location', 'northwest', 'Box', 'on', 'EdgeColor', 'k');
S.axes(ax); S.label(ax, 'A');

% ---- (B) and (C)
T = S.read('Fig5BC_titrations.csv');
lev = {'1x1x', '2x1x', '1x4x', '2x4x'}; lab = {'1\times, 1\times', '2\times, 1\times', '1\times, 4\times', '2\times, 4\times'};
cols = {S.k, S.r, S.b, S.g}; lw = [2.6 2.0 1.4 0.9];
ax = axes('Position', pos{2}); hold(ax, 'on');
for i = 1:4
    plot(ax, T.l, T.(['b_theta_' lev{i}]), '-', 'Color', cols{i}, 'LineWidth', lw(i), 'DisplayName', lab{i});
end
set(ax, 'XScale', 'log'); xlim(ax, [10^-1.2 10^1.2]); ylim(ax, [-0.02 1.02]);
xticks(ax, [0.1 1 10]); yticks(ax, [0 0.5 1]);
xlabel(ax, 'l'); ylabel(ax, '\theta'); legend(ax, 'Location', 'northeast', 'Box', 'on', 'EdgeColor', 'k');
S.axes(ax); S.label(ax, 'B');
ax = axes('Position', pos{3}); hold(ax, 'on');
for i = 1:4
    plot(ax, T.l, T.(['c_theta_bridge_' lev{i}]), '-', 'Color', cols{i}, 'LineWidth', lw(i));
end
set(ax, 'XScale', 'log'); xlim(ax, [10^-1.2 10^1.2]); ylim(ax, [-0.02 1.02]);
xticks(ax, [0.1 1 10]); yticks(ax, [0 0.5 1]);
xlabel(ax, 'l'); ylabel(ax, '\theta'); S.axes(ax); S.label(ax, 'C');

% ---- (D)
D = S.read('Fig5D_distributions.csv');
ax = axes('Position', pos{4}); hold(ax, 'on');
w = 0.38;
bar(ax, D.modified_sites-w/2, D.independent_nu1, w, 'FaceColor', S.b, 'EdgeColor', S.k, 'LineWidth', 0.6, 'DisplayName', '\nu = 1');
bar(ax, D.modified_sites+w/2, D.example_nu2p2, w, 'FaceColor', 'w', 'EdgeColor', S.k, 'LineWidth', 0.9, 'DisplayName', '\nu = 2.2');
xlim(ax, [-0.7 12.7]); ylim(ax, [0 0.25]); xticks(ax, [0 3 6 9 12]); yticks(ax, [0 0.1 0.2]);
xlabel(ax, 'i'); ylabel(ax, 'p_i'); legend(ax, 'Location', 'northeast', 'Box', 'on', 'EdgeColor', 'k');
S.axes(ax); set(ax, 'XMinorTick', 'off'); S.label(ax, 'D');
S.save(fig, 'Fig5');
