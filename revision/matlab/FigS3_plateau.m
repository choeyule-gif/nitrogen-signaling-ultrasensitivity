% S3 Fig.  Local coefficient of the total modified fraction as a surface over both totals
% (n = 3, A*Khat = 1, Bhat = Ghat = 0.4/Khat, K_m = 1/(n Bhat)); contours at 0.80, 0.90,
% 0.95 and 0.99 projected onto the floor.  Rotate interactively (rotate3d) and re-export.
S = plosstyle();
T = S.read('FigS3_plateau_surface.csv');
X = unique(T.TE_over_TS); Y = unique(T.TS_over_Km); nx = numel(X); ny = numel(Y);
Z = reshape(T.nH_local, nx, ny)';                 % rows: Y (T_S/K_m), columns: X (T_E/T_S)
lx = log10(X); ly = log10(Y);
fig = figure('Units', 'centimeters', 'Position', [2 2 11 8.5], 'Color', 'w');
ax = axes('Position', [0.02 0.05 0.96 0.92]); hold(ax, 'on');
surf(ax, lx, ly, Z, 'EdgeColor', [0.35 0.35 0.35], 'LineWidth', 0.25, 'FaceAlpha', 1);
colormap(ax, parula); clim(ax, [0.6 1.0]);
zfloor = 0.6;                                     % projected contours, drawn from contourc so that any release works
C = contourc(lx, ly, Z, [0.80 0.90 0.95 0.99]); k = 1;
while k < size(C, 2)
    m = C(2, k); seg = C(:, k+1:k+m);
    plot3(ax, seg(1, :), seg(2, :), zfloor*ones(1, m), 'k-', 'LineWidth', 0.7); k = k+m+1;
end
xlabel(ax, 'log_{10}(T_E/T_S)'); ylabel(ax, 'log_{10}(T_S/K_m)'); zlabel(ax, 'n_H^{loc}');
xticks(ax, [-4 -2 0]); yticks(ax, [-2 0 2]); zticks(ax, [0.6 0.8 1.0]); zlim(ax, [0.6 1.02]);
set(ax, 'FontSize', 8, 'FontName', 'Helvetica', 'LineWidth', 0.8, 'Box', 'off');
view(ax, -37.5, 18); grid(ax, 'off');
S.save(fig, 'FigS3');
