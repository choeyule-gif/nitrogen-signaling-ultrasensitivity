% Draws Fig1-Fig5 from ../figdata/*.csv into ../figures/ (vector PDF + PNG preview).
% Usage (MATLAB R2020a or later):   cd matlab; run_all_figures
here = fileparts(mfilename('fullpath')); addpath(here);
for k = 1:5
    fprintf('Fig%d ... ', k); run(fullfile(here, sprintf('Fig%d.m', k))); close all; fprintf('done\n');
end
