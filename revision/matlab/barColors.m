function colors = barColors()
%BARCOLORS  The two fill colours of the bar panels, taken from parula.
%   colors = barColors() returns a 2-by-3 matrix: a teal-leaning blue and an
%   orange-leaning yellow, read off the parula ramp at roughly 0.42 and 0.85.
%   Bars are large filled areas rather than thin lines, so they carry the
%   stronger colours of the S2 and S3 Fig surfaces directly and tie the bar
%   panels to those surfaces.  They are drawn without an outline.
%
%   See also SERIESCOLORS, PARULA.

colors = [
    0.13, 0.66, 0.86      % teal-leaning blue
    0.98, 0.75, 0.23];    % orange-leaning yellow
end
