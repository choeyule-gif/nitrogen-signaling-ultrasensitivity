function data = readFigData(fileName)
%READFIGDATA  Read one table of plotted values from ../figdata.
%   data = readFigData(fileName) returns ../figdata/<fileName> (relative to this
%   file) as a table.  Lines starting with # are comments and column names are
%   kept as written, for example 'Y_norm_PII_0.5'.  The files are produced by
%   ../export_figure_data.py and ../export_extra_figure_data.py; nothing is
%   recomputed while plotting.
%
%   See also READTABLE.

arguments
    fileName (1, :) char
end

dataDir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'figdata');
data = readtable(fullfile(dataDir, fileName), 'CommentStyle', '#', ...
    'VariableNamingRule', 'preserve');
end
