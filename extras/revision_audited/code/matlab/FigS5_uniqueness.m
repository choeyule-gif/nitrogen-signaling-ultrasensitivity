function FigS5_uniqueness()
%FIGS5_UNIQUENESS  S5 Fig: the residual whose zeros are the positive steady states.
%   FigS5A  both conditions of the uniqueness criterion satisfied: one zero,
%           R increasing
%   FigS5B  condition (i) broken by a spread of 2.89 in the catalytic ratios at
%           n = 4: three zeros
%   The zeros (sign changes of R) are marked with filled circles.  In FigS5B the
%   three zeros lie in a shallow fold of the residual, so that panel is drawn on
%   a tighter window, R from -0.3 to 0.35, which spreads the fold over a third
%   of the panel while still showing the residual rising away on both sides.
%
%   Reads figdata/FigS5A_residual.csv and figdata/FigS5B_residual.csv and
%   writes figures/FigS5A and figures/FigS5B (PDF + PNG).
%
%   See also RUN_ALL_FIGURES, NEWPANEL, SAVEPANEL, SERIESCOLORS.

%% FigS5A: criterion satisfied
residualA = readFigData('FigS5A_residual.csv');

[fig, ax] = newPanel();
plotResidual(ax, residualA);
xlim(ax, [1e-2, 1e2]);
ylim(ax, [-1.5, 6]);
savePanel(fig, 'FigS5A');

%% FigS5B: condition (i) broken, three zeros in a shallow fold
residualB = readFigData('FigS5B_residual.csv');

[fig, ax] = newPanel();
plotResidual(ax, residualB);
xlim(ax, [3e-4, 0.5]);
ylim(ax, [-0.3, 0.35]);
xticks(ax, [1e-3, 1e-2, 1e-1]);
savePanel(fig, 'FigS5B');
end

function plotResidual(ax, residual)
%PLOTRESIDUAL  Draw R(l) on a log axis with the zero line and the sign changes of R.
colors = seriesColors(2);
plot(ax, residual.l, residual.R, '-', 'Color', colors(1, :));
zeroCrossings = find(diff(sign(residual.R)) ~= 0);
plot(ax, residual.l(zeroCrossings), zeros(size(zeroCrossings)), 'o', ...
    'Color', colors(2, :), 'MarkerFaceColor', colors(2, :), 'MarkerSize', 6, ...
    'LineWidth', 0.75);
yline(ax, 0, 'k--', 'LineWidth', 1);
ax.XScale = 'log';
xlabel(ax, 'l');
ylabel(ax, 'R(l)');
end
