function Fig3()
%FIG3  Fig 3: the titration fit, one panel per file.
%   Fig3A  digitised titration (circles), the three bounded fixed-exponent
%          fits, and the pointwise 95% band of the free-exponent fit under the
%          coordinate perturbations
%   Fig3B  standardised residuals of the same three fits, with the band
%          |r| <= 2 sigma shaded
%   Fig3C  fitted exponent over the 2,000 perturbations, with the point fit
%          (solid) and the central 95% range (dotted)
%   Fig3D  profile likelihood loss against the regulatory swing, with the
%          criterion (dashed) and the lower endpoint (dotted)
%
%   The full manuscript also carries a logit panel and a likelihood-ratio panel
%   in its Fig 4; the shortened revision drops both, so neither is drawn here.
%
%   Reads figdata/Fig3*.csv and writes figures/Fig3A ... Fig3D (PDF + PNG).
%
%   See also RUN_ALL_FIGURES, NEWPANEL, SAVEPANEL, SERIESCOLORS.

exponents = [1, 2, 3];
exponentColors = seriesColors(numel(exponents));
baseColors = seriesColors(2);
blue = baseColors(1, :);
rose = baseColors(2, :);
grey = [0.5, 0.5, 0.5];
thinEdge = 0.75;                        % marker outlines, as in the reference figure

%% Fig3A: data, the three fixed-exponent fits, and the 95% band
points = readFigData('Fig3A_points.csv');
freeFit = readFigData('Fig3A_fit.csv');
fixedFits = readFigData('Fig3A_fixed_h_fits.csv');

[fig, ax] = newPanel();
fill(ax, [freeFit.glutamine_mM; flipud(freeFit.glutamine_mM)], ...
    [freeFit.band_lo; flipud(freeFit.band_hi)], paleTint(blue, 0.32), ...
    'EdgeColor', 'none', 'HandleVisibility', 'off');
for k = 1:numel(exponents)
    plot(ax, fixedFits.glutamine_mM, fixedFits.(sprintf('fit_h%d', exponents(k))), '-', ...
        'Color', exponentColors(k, :), 'DisplayName', sprintf('h = %d', exponents(k)));
end
plot(ax, points.glutamine_mM, points.UMP_per_trimer, 'ko', 'MarkerFaceColor', 'w', ...
    'LineWidth', thinEdge, 'HandleVisibility', 'off');
ax.XScale = 'log';
xlim(ax, [1e-2, 10^1.25]);
ylim(ax, [0, 3.15]);
xticks(ax, [1e-2, 1e-1, 1, 10]);
xlabel(ax, '[Gln] / mM');
ylabel(ax, 'UMP per trimer');
legend(ax, 'Location', 'southwest');
savePanel(fig, 'Fig3A');

%% Fig3B: standardised residuals of the same three fits
residuals = readFigData('Fig3B_residuals.csv');
glutamineLimits = [1e-2, 20];

[fig, ax] = newPanel();
fill(ax, glutamineLimits([1, 2, 2, 1]), [-2, -2, 2, 2], [0.9, 0.9, 0.9], ...
    'EdgeColor', 'none', 'HandleVisibility', 'off');
for k = 1:numel(exponents)
    plot(ax, residuals.glutamine_mM, residuals.(sprintf('r_h%d', exponents(k))), '-o', ...
        'Color', exponentColors(k, :), 'MarkerFaceColor', 'w', 'MarkerSize', 5, ...
        'LineWidth', thinEdge, 'DisplayName', sprintf('h = %d', exponents(k)));
end
yline(ax, 0, 'k-', 'HandleVisibility', 'off');
ax.XScale = 'log';
xlim(ax, glutamineLimits);
ylim(ax, [-8, 16]);                     % headroom so the legend clears the h = 1 peak
xticks(ax, [1e-2, 1e-1, 1, 10]);
yticks(ax, -8:4:16);
xlabel(ax, '[Gln] / mM');
ylabel(ax, 'r / \sigma');
legend(ax, 'Location', 'northeast', 'NumColumns', 3);
savePanel(fig, 'Fig3B');

%% Fig3C: fitted exponent over the perturbations
draws = readFigData('Fig3C_exponent_draws.csv');    % single column: the fitted h of each draw
markers = readFigData('Fig3C_markers.csv');
fittedExponent = markers.value(strcmp(markers.quantity, 'h_fit'));
lowerBound = markers.value(strcmp(markers.quantity, 'p025'));
upperBound = markers.value(strcmp(markers.quantity, 'p975'));

[fig, ax] = newPanel();
histogram(ax, draws{:, 1}, 35, 'Normalization', 'pdf', 'FaceColor', paleTint(blue, 0.55), ...
    'FaceAlpha', 1, 'EdgeColor', 'w');
xline(ax, fittedExponent, '-', 'Color', rose, 'LineWidth', 1.5);
xline(ax, lowerBound, ':', 'Color', grey, 'LineWidth', 1);
xline(ax, upperBound, ':', 'Color', grey, 'LineWidth', 1);
xlim(ax, [1.5, 3]);
ylim(ax, [0, 2.8]);
xlabel(ax, 'h');
ylabel(ax, 'density');
savePanel(fig, 'Fig3C');

%% Fig3D: profile likelihood loss against the regulatory swing
profileLoss = readFigData('Fig3D_swing_profile.csv');
markers = readFigData('Fig3D_markers.csv');
threshold = markers.value(strcmp(markers.quantity, 'threshold'));
lowerEndpoint = markers.value(strcmp(markers.quantity, 'lower95'));

[fig, ax] = newPanel();
plot(ax, profileLoss.swing, profileLoss.loss, '-', 'Color', blue);
yline(ax, threshold, '--', 'Color', grey, 'LineWidth', 1);
xline(ax, lowerEndpoint, ':', 'Color', rose, 'LineWidth', 1);
ax.XScale = 'log';
xlim(ax, [10^1.5, 1e5]);
ylim(ax, [0, 12]);
xticks(ax, [1e2, 1e3, 1e4, 1e5]);
xlabel(ax, 'R_{fwd} R_{rev}');
ylabel(ax, '2\Delta\itL');
savePanel(fig, 'Fig3D');
end
