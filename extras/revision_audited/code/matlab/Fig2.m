function Fig2()
%FIG2  Fig 2: the first layer, one panel per file.
%   Fig2A  the (nu, eta) plane.  Grey, the iso-coefficient hyperbolae
%          nu*eta = const, labelled by their value; blue, the hyperbola of the
%          comparison value 2.018 (not a measured local coefficient); orange, the ceiling that eq (6)
%          puts on eta, with the band above it shaded as unreachable by that
%          route; purple, the ceiling that eq (11) puts on nu at n = 3
%   Fig2B  maximum switch elasticity of eq (6) against the affinity ratio, with
%          the Mg (circle) and Mn (square) ratios
%   Fig2C  local logit slope of the fitted first-layer curve against glutamine,
%          the fitted exponent (dashed) and the value at the midpoint (circle)
%
%   Reads figdata/Fig2*.csv and writes figures/Fig2A, Fig2B, Fig2C (PDF + PNG).
%
%   See also RUN_ALL_FIGURES, NEWPANEL, SAVEPANEL, SERIESCOLORS.

colors = seriesColors(5);
blue = colors(1, :);
orange = colors(4, :);
purple = colors(5, :);
grey = [0.5, 0.5, 0.5];
lightGrey = [0.72, 0.72, 0.72];
thinEdge = 0.75;                        % marker outlines, as in the reference figure
etaCeiling = 2;

%% Fig2A: iso-coefficient curves, reference hyperbola and the two ceilings
isoCurves = readFigData('Fig2A_isocurves.csv');
hyperbola = readFigData('Fig2A_measured_hyperbola.csv');
nuLimits = [0.5, 3.1];
etaLimits = [0.5, 3.5];
isoValues = [1, 3, 4];                  % nu*eta = 2 is hidden under the reference 2.02 curve,
                                        % and nu*eta = 6 lies almost wholly in the shaded band
labelNu = 2.7;                          % iso-curve labels sit at this nu where they fit

[fig, ax] = newPanel();
% band above the eta ceiling, unreachable by the specified switch
fill(ax, nuLimits([1, 2, 2, 1]), [etaCeiling, etaCeiling, etaLimits(2), etaLimits(2)], ...
    paleTint(orange, 0.12), 'EdgeColor', 'none');
% curves of constant nu*eta, each labelled just above the curve with its value
for value = isoValues
    eta = isoCurves.(sprintf('eta_c%d', value));
    shown = eta >= etaLimits(1) & eta <= etaLimits(2);
    plot(ax, isoCurves.nu(shown), eta(shown), '-', 'Color', lightGrey, 'LineWidth', 0.75);
    nuHere = min(labelNu, value / (etaLimits(1) + 0.12));
    text(ax, nuHere, value / nuHere, num2str(value), 'FontSize', 8, 'Color', grey, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Margin', 1);
end
shown = hyperbola.eta >= etaLimits(1) & hyperbola.eta <= etaLimits(2);
plot(ax, hyperbola.nu(shown), hyperbola.eta(shown), '-', 'Color', blue);
yline(ax, etaCeiling, '-', 'Color', orange, 'LineWidth', 1.5);    % eta ceiling, eq (6)
xline(ax, 3, '-', 'Color', purple, 'LineWidth', 1.5);             % nu ceiling at n = 3, eq (11)
xlim(ax, nuLimits);
ylim(ax, etaLimits);
xticks(ax, [1, 2, 3]);
yticks(ax, [1, 2, 3]);
xlabel(ax, '\nu');
ylabel(ax, '\eta');
savePanel(fig, 'Fig2A');

%% Fig2B: maximum switch elasticity against the affinity ratio
etaMax = readFigData('Fig2B_eta_max.csv');
markers = readFigData('Fig2B_markers.csv');
markerShapes = {'o', 's'};    % Mg circle, Mn square

[fig, ax] = newPanel();
plot(ax, etaMax.affinity_ratio, etaMax.eta_max, '-', 'Color', orange, ...
    'HandleVisibility', 'off');
for k = 1:height(markers)
    plot(ax, markers.affinity_ratio(k), markers.eta_max(k), markerShapes{k}, ...
        'Color', 'k', 'MarkerFaceColor', 'w', 'LineWidth', thinEdge, ...
        'DisplayName', markers.ion{k});
end
yline(ax, etaCeiling, ':', 'Color', grey, 'LineWidth', 1, 'HandleVisibility', 'off');
ax.XScale = 'log';
xlim(ax, [1, 1e4]);
ylim(ax, [0.98, 2.04]);
xlabel(ax, 'A_{UT}/A_{UR}');
ylabel(ax, 'max \eta');
legend(ax, 'Location', 'southeast');
savePanel(fig, 'Fig2B');

%% Fig2C: local logit slope of the fitted first-layer curve
slope = readFigData('Fig2C_local_slope.csv');
markers = readFigData('Fig2C_markers.csv');
midpoint = markers.value(strcmp(markers.quantity, 'S_mM'));
slopeAtMidpoint = markers.value(strcmp(markers.quantity, 'nH_loc_at_S'));
fittedExponent = markers.value(strcmp(markers.quantity, 'nH_fit'));

[fig, ax] = newPanel();
plot(ax, slope.glutamine_mM, slope.nH_loc, '-', 'Color', blue);
plot(ax, midpoint, slopeAtMidpoint, 'o', 'Color', blue, 'MarkerFaceColor', 'w', ...
    'LineWidth', thinEdge);
yline(ax, fittedExponent, '--', 'Color', grey, 'LineWidth', 1);
ax.XScale = 'log';
xlim(ax, [1e-2, 10^1.3]);
ylim(ax, [0, 2.3]);
xticks(ax, [1e-2, 1e-1, 1, 10]);
xlabel(ax, '[Gln] / mM');
ylabel(ax, 'n_H');
savePanel(fig, 'Fig2C');
end
