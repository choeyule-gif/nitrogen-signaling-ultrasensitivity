function S = plosstyle()
%PLOSSTYLE  Shared settings for the manuscript figures (MATLAB look, PLOS sizes).
%   S = plosstyle() returns a struct with colours, the figdata folder, the output
%   folder and helper function handles:
%     S.read(name)          table read from figdata/<name> (comment lines skipped)
%     S.axes(ax)            box on, ticks in, minor ticks, line widths, font size
%     S.label(ax, 'A')      bold panel letter outside the top-left corner
%     S.save(fig, 'Fig1')   vector PDF (and PNG preview) in the figures folder
%   Nothing is recomputed here: every plotted number comes from a CSV written by
%   ../export_figure_data.py.
here = fileparts(mfilename('fullpath'));
S.figdata = fullfile(here, '..', 'figdata');
S.outdir  = fullfile(here, '..', 'figures');
if ~exist(S.outdir, 'dir'), mkdir(S.outdir); end
% colours (black, red, blue, green, grey, orange, purple)
S.k  = [0 0 0];        S.r = [0.91 0 0.04];   S.b = [0 0 0.93];
S.g  = [0 0.63 0];     S.gr = [0.5 0.5 0.5];  S.o = [0.85 0.33 0.10];
S.p  = [0.49 0.18 0.56]; S.lb = [0.80 0.87 0.96];
S.fs = 8;                                   % tick label font size
S.read  = @(name) readtable(fullfile(S.figdata, name), 'CommentStyle', '#', 'VariableNamingRule', 'preserve');
S.axes  = @style_axes;
S.label = @panel_label;
S.save  = @(fig, name) save_fig(fig, fullfile(S.outdir, name));
end

function style_axes(ax)
set(ax, 'Box', 'on', 'TickDir', 'in', 'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'LineWidth', 0.9, 'FontSize', 8, 'Layer', 'top', 'FontName', 'Helvetica');
ax.XLabel.FontSize = 9; ax.YLabel.FontSize = 9;
end

function panel_label(ax, letter)
text(ax, -0.24, 1.07, letter, 'Units', 'normalized', 'FontWeight', 'bold', ...
    'FontSize', 11, 'FontName', 'Helvetica', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
end

function save_fig(fig, path)
exportgraphics(fig, [path '.pdf'], 'ContentType', 'vector');
exportgraphics(fig, [path '.png'], 'Resolution', 220);
end
