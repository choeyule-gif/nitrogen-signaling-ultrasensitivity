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

outDir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'figures', 'matlab');
if ~isfolder(outDir)
    mkdir(outDir);
end
S=zstyle();
axs=findall(fig,'Type','axes');
for a=reshape(axs,1,[])
 zax(a,S); colormap(a,S.seqBlue);
end
for l=reshape(findall(fig,'Type','legend'),1,[])
 set(l,'FontName',S.font,'FontSize',10,'Box','off','TextColor',S.ink1);
end
exportgraphics(fig, fullfile(outDir, [name, '.pdf']), 'ContentType', 'vector','BackgroundColor','white');
exportgraphics(fig, fullfile(outDir, [name, '.png']), 'Resolution', 600);
end
