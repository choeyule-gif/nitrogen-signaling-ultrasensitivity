% Draws Fig1-Fig5 from ../figdata/*.csv into ../figures/ (vector PDF + PNG preview).
% Usage (MATLAB R2020a or later):   cd matlab; run_all_figures
here = fileparts(mfilename('fullpath')); addpath(here);
for k = 1:5
    fprintf('Fig%d ... ', k); run(fullfile(here, sprintf('Fig%d.m', k))); close all; fprintf('done\n');
end
for f = {'FigS1_invariance', 'FigS2_nu_ceiling', 'FigS3_plateau', 'FigS4_nu2', 'FigS5_uniqueness'}
    fprintf('%s ... ', f{1}); run(fullfile(here, [f{1} '.m'])); close all; fprintf('done\n');
end
