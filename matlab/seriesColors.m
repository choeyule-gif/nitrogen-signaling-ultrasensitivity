function colors = seriesColors(nSeries)
%SERIESCOLORS  Muted categorical line colours, in the style of Nature figures.
%   colors = seriesColors(nSeries) returns the first nSeries rows of a fixed
%   seven-colour palette: blue, rose, teal, orange, purple, green, grey.
%
%   The palette is pastel but not pale.  Every entry sits between 0.48 and 0.65
%   relative luminance, so a 1.5 pt line and a 0.75 pt marker edge both hold up
%   on white, and no entry shouts over the others.  Blue, teal, orange and green
%   are drawn from the range the parula surfaces of S2 and S3 Figs cover, so the
%   line panels and the surface panels sit together; rose and purple are there
%   to separate series that parula-like hues alone would not.
%
%   The order puts the two highest-contrast entries first, so a two-series panel
%   gets blue against rose.
%
%   See also NEWPANEL, PALETINT.

arguments
    nSeries (1, 1) double {mustBePositive, mustBeInteger}
end

palette = [
    0.31, 0.51, 0.74     % blue
    0.84, 0.47, 0.45     % rose
    0.33, 0.66, 0.64     % teal
    0.90, 0.58, 0.22     % orange
    0.58, 0.50, 0.74     % purple
    0.52, 0.70, 0.40     % green
    0.55, 0.55, 0.55];   % grey

if nSeries > size(palette, 1)
    error('seriesColors:tooManySeries', ...
        'The palette holds %d colours; %d were asked for.', size(palette, 1), nSeries);
end
colors = palette(1:nSeries, :);
end
