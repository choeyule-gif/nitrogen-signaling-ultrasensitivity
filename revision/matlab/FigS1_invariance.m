function FigS1_invariance()
%FIGS1_INVARIANCE  S1 Fig: invariance of the ladder ratios, one panel per file.
%   FigS1A  ratios s_i/s_{i-1} at n = 4 (non-identical ladder) against the free
%           enzyme over a 5,000-fold range, solved in the raw species variables
%           without imposing eq (1)
%   FigS1B  the same ratios at n = 3 while a dead-end complex E.S_0 has its
%           association constant swept over eleven decades (T_E/T_L = 1.1e-9)
%
%   Every ratio is flat, which is the point of the figure.  The vertical range
%   carries a spare decade above the highest ratio so that the legend sits in
%   empty space instead of over the curves.
%
%   Reads figdata/FigS1A_enzyme_sweep.csv and figdata/FigS1B_deadend_sweep.csv
%   and writes figures/FigS1A and figures/FigS1B (PDF + PNG).
%
%   See also RUN_ALL_FIGURES, NEWPANEL, SAVEPANEL, SERIESCOLORS.

%% FigS1A: ratios against the free enzyme
enzymeSweep = readFigData('FigS1A_enzyme_sweep.csv');
colors = seriesColors(4);

[fig, ax] = newPanel();
for k = 1:4
    plot(ax, enzymeSweep.free_enzyme, enzymeSweep.(sprintf('s%d_over_s%d', k, k - 1)), ...
        '-', 'Color', colors(k, :), 'DisplayName', sprintf('s_%d/s_%d', k, k - 1));
end
ax.XScale = 'log';
ax.YScale = 'log';
xlim(ax, [min(enzymeSweep.free_enzyme), max(enzymeSweep.free_enzyme)]);
ylim(ax, [0.1, 300]);                  % headroom above the top ratio for the legend
yticks(ax, [0.1, 1, 10]);
xlabel(ax, 'e');
ylabel(ax, 's_i/s_{i-1}');
legend(ax, 'Location', 'north', 'NumColumns', 2);
savePanel(fig, 'FigS1A');

%% FigS1B: ratios against the dead-end association constant
deadEndSweep = readFigData('FigS1B_deadend_sweep.csv');
colors = seriesColors(3);

[fig, ax] = newPanel();
for k = 1:3
    plot(ax, deadEndSweep.K_Z, deadEndSweep.(sprintf('s%d_over_s%d', k, k - 1)), ...
        '-', 'Color', colors(k, :), 'DisplayName', sprintf('s_%d/s_%d', k, k - 1));
end
ax.XScale = 'log';
ax.YScale = 'log';
xlim(ax, [1e-6, 1e5]);
ylim(ax, [0.25, 60]);                  % headroom above the top ratio for the legend
yticks(ax, [0.3, 1, 3]);
xlabel(ax, 'K_Z');
ylabel(ax, 's_i/s_{i-1}');
legend(ax, 'Location', 'north', 'NumColumns', 3);
savePanel(fig, 'FigS1B');
end
