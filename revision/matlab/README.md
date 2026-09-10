# matlab/

MATLAB (R2021b or later: `exportgraphics` needs R2020a, `xline(ax, ...)` needs R2021b)
scripts for the figures of the revised manuscript and for the companion figure.  Every plotted value is read from `../figdata/*.csv`, written by
`../export_figure_data.py` and `../export_extra_figure_data.py`; nothing is recomputed
while plotting.

```matlab
cd matlab
run_all_figures        % Fig 1-5 and S1-S5 Figs, one file per panel, into ../figures/
Fig_ReF_companion      % the Re F(i omega) figure of the companion stability manuscript
```

## One panel, one file

Each panel is drawn in its own figure window and exported on its own, as a vector PDF
plus a 300-dpi PNG preview (`Fig2A.pdf`, `Fig2B.pdf`, `Fig2C.pdf`, ...).  The panels
carry no letters and no titles; letters are added when the panels are assembled.
Fig 1 no longer has the schematic panel: `Fig1A` is the theta(l) panel and `Fig1B` the
local-coefficient panel (they read `Fig1B_theta.csv` and `Fig1C_local_hill.csv`, whose
names keep the letters of the old layout).

| script | panels written |
| --- | --- |
| `Fig1.m` | `Fig1A`, `Fig1B` |
| `Fig2.m` | `Fig2A`, `Fig2B`, `Fig2C` |
| `Fig3.m` | `Fig3A`, `Fig3B`, `Fig3C`, `Fig3D` |
| `Fig4.m` | `Fig4A`, `Fig4B`, `Fig4C` |
| `Fig5.m` | `Fig5A`, `Fig5B`, `Fig5C`, `Fig5D` |
| `FigS1_invariance.m` | `FigS1A`, `FigS1B` |
| `FigS2_nu_ceiling.m` | `FigS2` |
| `FigS3_plateau.m` | `FigS3` |
| `FigS4_nu2.m` | `FigS4` |
| `FigS5_uniqueness.m` | `FigS5A`, `FigS5B` |
| `Fig_ReF_companion.m` | `Fig_ReF_companion` |

## Style

The frame follows the figures of Nature: only the left and bottom rules are drawn, ticks
point outwards, and legends carry no box.  The rest is the MATLAB default, Helvetica
10 pt with axis labels at 11 pt, 0.5 pt axis lines, minor ticks on logarithmic axes only,
no grid on 2-D panels (grid on for the surfaces, as `surf` draws it).  Panels are
8.8 cm x 6.6 cm, 11 cm x 8.5 cm for the surfaces, and data lines are 1.5 pt.

Series colours come from `seriesColors(n)`, a fixed seven-colour palette: blue, rose,
teal, orange, purple, green, grey.  It is pastel but not pale, every entry between 0.48 and
0.65 relative luminance, so thin lines and 0.75 pt marker edges hold up on white and no
entry shouts over the others.  Blue, teal, orange and green come from the range the parula
surfaces of S2 and S3 Figs cover, so the line panels and the surface panels sit together;
rose and purple separate series that parula-like hues alone would not.  The order puts
the two highest-contrast entries first, so a two-series panel gets blue against rose.

The surfaces of `FigS2` and `FigS3` keep `parula`: a sequential colormap suits a height
field.  Only the curves drawn over the S2 Fig surface take the shared palette.

The two bar panels, `Fig4B` and `Fig5D`, are filled from `barColors`, a teal-leaning blue
and an orange-leaning yellow read straight off the parula ramp, and are drawn without an
outline.  Bars are large filled areas rather than thin lines, so they carry the stronger
surface colours directly.

Curves that coincide also carry different line styles, since colour alone would hide all
but the last one drawn.  Reference lines are grey or black `xline`/`yline`s and shaded
regions are opaque tints from `paleTint`, so the PDFs contain no transparency and stay
pure vector.

## Choosing a viewpoint for S3 Fig

`FigS3_plateau` takes an optional N-by-2 matrix of `[azimuth, elevation]` pairs.  With
several rows it writes one sample per row, named `FigS3_view<k>_az<az>_el<el>`, so a
viewpoint can be picked by eye instead of by rotating:

```matlab
FigS3_plateau([-37.5 18; -135 25; -60 40; 30 20; -100 55; -20 12])
```

Called with no argument it writes `FigS3` at the default viewpoint, which is what
`run_all_figures` does.  `Fig_ReF_companion` takes the same argument, and writes its samples
as `Fig_ReF_az<az>_el<el>`.

Helpers: `newPanel` (figure and axes in the shared style), `readFigData` (one CSV table
from `../figdata`), `savePanel` (PDF and PNG into `../figures`), `seriesColors` (the
palette), `paleTint` (opaque fills).  Apart from `FigS3_plateau` each `Fig*.m` is a
function without arguments.  After rotating a surface interactively, re-export it with
`savePanel(gcf, 'FigS2')`.
