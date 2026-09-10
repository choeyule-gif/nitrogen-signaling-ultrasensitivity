function tinted = paleTint(color, strength)
%PALETINT  Mix a colour with white, for opaque fills behind data.
%   tinted = paleTint(color, strength) keeps the given fraction of the colour,
%   so strength = 0.2 is a light tint and strength = 1 is the colour itself.
%   The result is opaque, which keeps the exported PDF free of transparency.
%
%   See also SERIESCOLORS.

arguments
    color (1, 3) double
    strength (1, 1) double {mustBeNonnegative, mustBeLessThanOrEqual(strength, 1)}
end

tinted = 1 - strength * (1 - color);
end
