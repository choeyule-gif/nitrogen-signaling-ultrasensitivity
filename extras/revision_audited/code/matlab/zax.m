function zax(ax, S)
%ZAX  Journal-style axes (like the reference isotherm panel): no box, outward ticks,
%   recessive hairline axes in secondary ink, system sans, no grid.
if nargin < 2, S = zstyle(); end
set(ax, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.012 0.012], 'LineWidth', 0.9, ...
    'FontName', S.font, 'FontSize', S.fs, 'XColor', S.ink2, 'YColor', S.ink2, 'ZColor', S.ink2, ...
    'Color', 'none', 'XGrid', 'off', 'YGrid', 'off', 'Layer', 'top');
ax.XLabel.Color = S.ink1; ax.YLabel.Color = S.ink1; ax.ZLabel.Color = S.ink1;
ax.Title.Color = S.ink1; ax.Title.FontWeight = 'normal'; ax.Title.FontSize = S.fs + 1;
end
