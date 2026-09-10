function FigS2_nu_ceiling()
%FIGS2_NU_CEILING  S2 Fig: supremum of the ladder factor nu over three-site ladders.
%   Surface of sup_u nu for every three-site ladder c = (1, 3a, 3b, 1) over
%   (log10 a, log10 b).  a = b = 1 is the binomial ladder (nu = 1, circle) and
%   the ceiling nu = 3 is approached as a, b -> 0.  The two lines are the
%   ladders that reproduce the observed ratio 1.1106 with the fitted-exponent
%   (blue) and range (orange) estimators.  Rotate interactively, then re-export
%   with savePanel(gcf, 'FigS2').
%
%   Reads figdata/FigS2_*.csv and writes figures/FigS2 (PDF + PNG).
%
%   See also RUN_ALL_FIGURES, NEWPANEL, SAVEPANEL.

surfaceTable = readFigData('FigS2_nu_sup_surface.csv');
locusFit = readFigData('FigS2_locus_fit.csv');
locusRange = readFigData('FigS2_locus_range.csv');
log10a = unique(surfaceTable.log10_a);
log10b = unique(surfaceTable.log10_b);
nuSup = reshape(surfaceTable.nu_sup, numel(log10a), numel(log10b))';    % rows: log10 b

[fig, ax] = newPanel(11, 8.5);
colors = seriesColors(2);
surf(ax, log10a, log10b, nuSup, 'EdgeColor', [0.35, 0.35, 0.35], 'LineWidth', 0.25);
plot3(ax, 0, 0, 1, 'ko', 'MarkerFaceColor', 'w');                       % binomial ladder
plot3(ax, locusFit.log10_a, locusFit.log10_b, locusFit.nu_sup + 0.02, '-', ...
    'Color', colors(1, :));
plot3(ax, locusRange.log10_a, locusRange.log10_b, locusRange.nu_sup + 0.02, '-', ...
    'Color', colors(2, :));
ax.CLim = [1, 3];
zlim(ax, [1, 3.05]);
xlabel(ax, 'log_{10} a');
ylabel(ax, 'log_{10} b');
zlabel(ax, 'sup_u \nu');
grid(ax, 'on');
view(ax, -40, 24);
savePanel(fig, 'FigS2');
end
