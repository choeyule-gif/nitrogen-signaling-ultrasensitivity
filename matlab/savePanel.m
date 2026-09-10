function savePanel(fig, name)
%SAVEPANEL  Export one figure as a vector PDF and a 300 dpi PNG preview.
%   savePanel(fig, name) writes ../figures/<name>.pdf and ../figures/<name>.png
%   (relative to this file), creating the folder if needed.
%
%   See also NEWPANEL, EXPORTGRAPHICS.

arguments
    fig (1, 1) matlab.ui.Figure
    name (1, :) char
end

outDir = fullfile(fileparts(mfilename('fullpath')), '..', 'figures');
if ~isfolder(outDir)
    mkdir(outDir);
end
exportgraphics(fig, fullfile(outDir, [name, '.pdf']), 'ContentType', 'vector');
exportgraphics(fig, fullfile(outDir, [name, '.png']), 'Resolution', 300);
end
