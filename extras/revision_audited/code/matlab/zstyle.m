function S = zstyle()
%ZSTYLE  Shared style for every figure: fonts, ink, species palette, element colors.
%   Saturated author palette; retain line styles and labels for redundant encoding.
%   Palette sampled from the author Figure 1 after ICC-to-sRGB conversion.
%   Cu (amber, 2.1:1 on white) always carries a direct label or a table value (contrast relief).
S.font = 'Helvetica';
S.fs = 11;
S.ink1 = hex('#0b0b0b');      % primary ink
S.ink2 = hex('#52514e');      % secondary ink (axes, ticks)
S.ink3 = hex('#898781');      % muted ink
S.grid = hex('#e1e0d9');
S.species = {'Na', 'Mn', 'Co', 'Ni', 'Cu', 'Zn'};
S.col.Na = hex('#8a8984');
S.col.Mn = hex('#0690f6');
S.col.Co = hex('#ff7b02');
S.col.Ni = hex('#2bd302');
S.col.Cu = hex('#ff1f1f');
S.col.Zn = hex('#0e56ff');
S.marker.Na = 'o'; S.marker.Mn = 's'; S.marker.Co = 'd'; S.marker.Ni = '^'; S.marker.Cu = 'o'; S.marker.Zn = 'v';
S.label.Na = 'Na^+'; S.label.Mn = 'Mn^{2+}'; S.label.Co = 'Co^{2+}'; S.label.Ni = 'Ni^{2+}';
S.label.Cu = 'Cu^{2+}'; S.label.Zn = 'Zn^{2+}';
% three ordered conditions (reference isotherm style): blue / red / aqua, validated all-pairs
S.cond = [hex('#0690f6'); hex('#ff7b02'); hex('#2bd302')];
% framework atoms (structure renderings): O red, Si pale gold (kept apart from Cu amber), Al grey-lilac
S.el.O = hex('#e04b3c'); S.el.Si = hex('#e6d28c'); S.el.Al = hex('#b5a2c8');
S.rad.O = 0.30; S.rad.Si = 0.38; S.rad.Al = 0.40; S.rad.cation = 0.95; S.rad.Na = 0.80;
% sequential ramps (one hue, light -> dark) for densities
S.seqBlue = interp1([0 0.25 0.5 0.75 1], [hex('#fcfcfb'); hex('#d0e8fd'); hex('#7abff9'); hex('#0690f6'); hex('#054778')], linspace(0, 1, 256));
S.seqOrange = interp1([0 0.25 0.5 0.75 1], [hex('#fcfcfb'); hex('#fbd3bf'); hex('#f39a6d'); hex('#eb6834'); hex('#9c3a12')], linspace(0, 1, 256));
S.seqGray = interp1([0 0.5 1], [hex('#fcfcfb'); hex('#b9b8b2'); hex('#52514e')], linspace(0, 1, 256));
end

function c = hex(h)
h = strrep(h, '#', '');
c = [hex2dec(h(1:2)) hex2dec(h(3:4)) hex2dec(h(5:6))] / 255;
end
