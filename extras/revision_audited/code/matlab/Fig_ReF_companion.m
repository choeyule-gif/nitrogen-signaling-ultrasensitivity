function Fig_ReF_companion(viewAngles)
%FIG_REF_COMPANION  Companion (stability) manuscript: Re F(i omega) across n.
%   For each n the curve Re F(i omega) is drawn at its own depth in the
%   (log omega, n, Re F) box, and the two frequency ranges are marked as filled
%   regions along the floor of the box rather than as lengths of line.  Blue is
%   the range the proofs cover: every omega at n = 1, and omega above the
%   circulation threshold otherwise.  Rose is the band that remains open,
%   n >= 2 below the threshold.  Circles at omega -> 0 mark F(0).
%
%   The regions are low ribbons, not full curtains under the curves.  Curtains
%   were the first thing tried and do not work here: they reach almost the top
%   of the box, they have to stay opaque to keep the PDF vector, and at every
%   viewpoint the nearer ones bury the farther ones, so five of the six ranges
%   disappear.  A ribbon shows the same interval and hides nothing.
%
%   Fig_ReF_companion(viewAngles) takes an N-by-2 matrix of
%   [azimuth, elevation] pairs.  With one row it writes figures/Fig_ReF_companion;
%   with several it writes one sample per row, named Fig_ReF_az<az>_el<el>, so a
%   viewpoint can be picked by eye.
%
%   Reads figdata/Fig_ReF_companion.csv and figdata/Fig_ReF_markers.csv.
%
%   See also NEWPANEL, SAVEPANEL, SERIESCOLORS, PALETINT.

arguments
    viewAngles (:, 2) double = [-50, 30]
end

ribbonHeight = 0.12;                    % height of the floor ribbons, in Re F
response = readFigData('Fig_ReF_companion.csv');
markers = readFigData('Fig_ReF_markers.csv');
siteCounts = unique(response.n);
colors = seriesColors(2);
proved = colors(1, :);
open = colors(2, :);
provedFill = paleTint(proved, 0.30);
openFill = paleTint(open, 0.38);
grey = [0.5, 0.5, 0.5];

[fig, ax] = newPanel(11, 8.5);
for k = numel(siteCounts):-1:1          % back of the box first
    n = siteCounts(k);
    rows = response.n == n;
    logOmega = log10(response.omega(rows));
    reF = response.ReF(rows);
    firstProved = find(response.omega(rows) >= markers.circulation_threshold(markers.n == n), 1);
    below = 1:firstProved;
    above = firstProved:numel(logOmega);
    if n == 1                           % proved for every omega
        fillRibbon(ax, logOmega, n, ribbonHeight, provedFill);
    else
        fillRibbon(ax, logOmega(below), n, ribbonHeight, openFill);
        fillRibbon(ax, logOmega(above), n, ribbonHeight, provedFill);
        plot3(ax, logOmega(below), n * ones(size(below)), reF(below), '-', ...
            'Color', open, 'LineWidth', 2.4);
    end
    plot3(ax, logOmega(above), n * ones(size(above)), reF(above), '-', 'Color', proved);
    if n == 1
        plot3(ax, logOmega(below), n * ones(size(below)), reF(below), '-', 'Color', proved);
    end
    plot3(ax, logOmega(1), n, markers.F0(markers.n == n), 'o', 'Color', proved, ...
        'MarkerFaceColor', 'w', 'MarkerSize', 4, 'LineWidth', 0.75);
end
plot3(ax, [-3, 3], [0.3, 0.3], [0, 0], '-', 'Color', grey, 'LineWidth', 0.75);
xticks(ax, [-3, -1, 1, 3]);
yticks(ax, siteCounts);
zticks(ax, [0, 0.5, 1]);
zlim(ax, [0, 1.05]);
xlabel(ax, 'log_{10}\omega');
ylabel(ax, 'n');
zlabel(ax, 'Re F');
grid(ax, 'on');

for k = 1:size(viewAngles, 1)
    view(ax, viewAngles(k, 1), viewAngles(k, 2));
    if size(viewAngles, 1) == 1
        savePanel(fig, 'Fig_ReF_companion');
    else
        savePanel(fig, sprintf('Fig_ReF_az%+d_el%+d', round(viewAngles(k, 1)), ...
            round(viewAngles(k, 2))));
    end
end
end

function fillRibbon(ax, logOmega, depth, height, faceColor)
%FILLRIBBON  Fill a ribbon of the given height along the floor, at one depth in n.
logOmega = logOmega(:);
edges = [min(logOmega); max(logOmega)];
fill3(ax, edges([1, 2, 2, 1]), depth * ones(4, 1), [0; 0; height; height], ...
    faceColor, 'EdgeColor', 'none');
end
