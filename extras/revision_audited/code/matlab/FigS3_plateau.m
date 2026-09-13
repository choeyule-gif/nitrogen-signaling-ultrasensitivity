function FigS3_plateau(viewAngles)
%FIGS3_PLATEAU  S3 Fig: local coefficient of the total modified fraction.
%   FigS3_plateau() draws the surface of the local coefficient over
%   log10(T_E/T_S) and log10(T_S/K_m) (n = 3, A*Khat = 1, Bhat = Ghat =
%   0.4/Khat, K_m = 1/(n Bhat)) at the default viewpoint and writes
%   figures/FigS3.  The contours at 0.80, 0.90, 0.95 and 0.99 are projected onto
%   the floor.
%
%   FigS3_plateau(viewAngles) takes an N-by-2 matrix of [azimuth, elevation]
%   pairs in degrees.  With one row it writes figures/FigS3 as above; with
%   several it writes one sample per row, named FigS3_view<k>_az<az>_el<el>, so
%   that a viewpoint can be chosen by eye.  For example
%
%       FigS3_plateau([-37.5 18; -135 25; -60 40; 30 20])
%
%   Reads figdata/FigS3_plateau_surface.csv.
%
%   See also RUN_ALL_FIGURES, NEWPANEL, SAVEPANEL, VIEW.

arguments
    viewAngles (:, 2) double = [-37.5, 18]
end

floorLevel = 0.6;
contourLevels = [0.80, 0.90, 0.95, 0.99];

surfaceTable = readFigData('FigS3_plateau_surface.csv');
enzymeRatio = log10(unique(surfaceTable.TE_over_TS));
substrateRatio = log10(unique(surfaceTable.TS_over_Km));
localHill = reshape(surfaceTable.nH_local, numel(enzymeRatio), numel(substrateRatio))';

[fig, ax] = newPanel(11, 8.5);
surf(ax, enzymeRatio, substrateRatio, localHill, 'EdgeColor', [0.35, 0.35, 0.35], ...
    'LineWidth', 0.25);
% contours projected onto the floor, taken from contourc so that any release works
contourMatrix = contourc(enzymeRatio, substrateRatio, localHill, contourLevels);
k = 1;
while k < size(contourMatrix, 2)
    nPoints = contourMatrix(2, k);
    segment = contourMatrix(:, k + 1:k + nPoints);
    plot3(ax, segment(1, :), segment(2, :), floorLevel * ones(1, nPoints), 'k-', ...
        'LineWidth', 0.75);
    k = k + nPoints + 1;
end
ax.CLim = [floorLevel, 1];
zlim(ax, [floorLevel, 1.02]);
xlabel(ax, 'log_{10}(T_E/T_S)');
ylabel(ax, 'log_{10}(T_S/K_m)');
zlabel(ax, 'n_H^{loc}');
grid(ax, 'on');

for k = 1:size(viewAngles, 1)
    view(ax, viewAngles(k, 1), viewAngles(k, 2));
    if size(viewAngles, 1) == 1
        savePanel(fig, 'FigS3');
    else
        savePanel(fig, sprintf('FigS3_view%d_az%+d_el%+d', k, ...
            round(viewAngles(k, 1)), round(viewAngles(k, 2))));
    end
end
end
