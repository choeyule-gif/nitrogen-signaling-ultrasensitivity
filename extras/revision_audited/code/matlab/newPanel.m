function [fig, ax] = newPanel(widthCm, heightCm)
%NEWPANEL  Create a blank single-panel figure in the shared manuscript style.
%   [fig, ax] = newPanel() returns an 8.8 cm x 6.6 cm white figure holding one
%   axes.  [fig, ax] = newPanel(widthCm, heightCm) sets the size in centimetres.
%
%   The frame follows the figures of Nature: only the left and bottom rules are
%   drawn, ticks point outwards, and legends carry no box.  The rest is the
%   MATLAB default (Helvetica 10 pt, 0.5 pt axis lines, minor ticks on
%   logarithmic axes only), with hold on so that a panel is built up by
%   successive plot calls and 1.5 pt data lines.
%
%   From R2025a the light theme is pinned.  Without that the figure follows the
%   desktop appearance, so the same script exports white panels by day and dark
%   panels at night; on earlier releases every figure is light already.
%
%   See also SAVEPANEL, READFIGDATA, SERIESCOLORS.

arguments
    widthCm (1, 1) double {mustBePositive} = 11
    heightCm (1, 1) double {mustBePositive} = 10
end

fig = figure('Units', 'centimeters', 'Position', [2, 2, widthCm, heightCm], ...
    'Color', 'w', 'DefaultLineLineWidth', 1.5, 'DefaultLegendBox', 'off');
if isprop(fig, 'Theme')
    fig.Theme = 'light';
end
ax = axes(fig);
hold(ax, 'on');
box(ax, 'off');            % no top or right rule
ax.TickDir = 'out';
S=zstyle();zax(ax,S);set(fig,'Visible','off','DefaultLineLineWidth',2);

end
