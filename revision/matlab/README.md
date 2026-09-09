# matlab/

MATLAB (R2020a or later; uses `exportgraphics`) scripts for the five figures of the revised
manuscript.  Every plotted value is read from `../figdata/*.csv`, written by
`../export_figure_data.py`; nothing is recomputed while plotting.

```matlab
cd matlab
run_all_figures          % writes ../figures/Fig1.pdf ... Fig5.pdf and FigS1 ... FigS5 (+ PNG previews)
Fig_ReF_companion        % the Re F(i omega) figure of the companion stability manuscript
```

`plosstyle.m` holds the shared style (box on, ticks inward, minor ticks, 8-pt tick labels,
black/red/blue/green/grey palette, open markers) and the panel-letter and export helpers.
The figures carry no prose: numbers, symbols and legends only.
