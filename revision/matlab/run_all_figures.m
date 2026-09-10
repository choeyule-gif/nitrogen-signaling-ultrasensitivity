%RUN_ALL_FIGURES  Draw every panel of Fig 1-5 and S1-S5 Figs from ../figdata.
%   Writes one vector PDF and one 300 dpi PNG per panel into ../figures.
%   Usage (MATLAB R2020a or later):   cd matlab; run_all_figures
%   The companion figure is drawn on its own:   Fig_ReF_companion

here = fileparts(mfilename('fullpath'));
addpath(here);
figureFunctions = {@Fig1, @Fig2, @Fig3, @Fig4, @Fig5, @FigS1_invariance, ...
    @FigS2_nu_ceiling, @FigS3_plateau, @FigS4_nu2, @FigS5_uniqueness};
for k = 1:numel(figureFunctions)
    fprintf('%s ... ', func2str(figureFunctions{k}));
    figureFunctions{k}();
    close all
    fprintf('done\n');
end
