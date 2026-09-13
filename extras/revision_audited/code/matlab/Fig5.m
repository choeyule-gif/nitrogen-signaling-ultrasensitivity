function Fig5()
%FIG5  Fig 5: binding laws, robustness of the titration and site distributions.
%   Fig5A  ligands bound per enzyme for a one-site (h = 1) and a concerted
%          two-site (h = 2) binding law
%   Fig5B  glutamine titrations predicted by the minimal network at 1x/2x
%          enzyme and 1x/4x target; the four curves coincide
%   Fig5C  the same conditions with a catalytically active complex bridging
%          two targets
%   Fig5D  twelve-site distributions with the same mean for nu = 1 and nu = 2.2
%
%   Reads figdata/Fig5*.csv and writes figures/Fig5A ... Fig5D (PDF + PNG).
%
%   See also RUN_ALL_FIGURES, NEWPANEL, SAVEPANEL, SERIESCOLORS.

%% Fig5A: one-site and concerted two-site binding
binding = readFigData('Fig5A_binding.csv');
colors = seriesColors(2);

[fig, ax] = newPanel();
plot(ax, binding.scaled_glutamine, binding.one_ligand, '-', 'Color', colors(1, :), ...
    'DisplayName', 'h = 1');
plot(ax, binding.scaled_glutamine, binding.two_ligands, '-', 'Color', colors(2, :), ...
    'DisplayName', 'h = 2');
ax.XScale = 'log';
xlim(ax, [1e-2, 1e2]);
ylim(ax, [0, 2.05]);
xlabel(ax, '[Gln] / K_L');
ylabel(ax, 'bound / enzyme');
legend(ax, 'Location', 'northwest');
savePanel(fig, 'Fig5A');

%% Fig5B and Fig5C: titrations at 1x/2x enzyme and 1x/4x target
titrations = readFigData('Fig5BC_titrations.csv');

[fig, ax] = newPanel();
plotTitrations(ax, titrations, 'b_theta_');            % minimal network
savePanel(fig, 'Fig5B');

[fig, ax] = newPanel();
plotTitrations(ax, titrations, 'c_theta_bridge_');     % bridging complex
savePanel(fig, 'Fig5C');

%% Fig5D: twelve-site distributions with the same mean
distributions = readFigData('Fig5D_distributions.csv');
colors = barColors();

[fig, ax] = newPanel();
plot(ax,distributions.modified_sites,distributions.independent_nu1,'-o','Color',colors(1,:),'MarkerFaceColor',colors(1,:),'LineWidth',2,'DisplayName','\nu = 1');
plot(ax,distributions.modified_sites,distributions.example_nu2p2,'-s','Color',colors(2,:),'MarkerFaceColor',colors(2,:),'LineWidth',2,'DisplayName','\nu = 2.2');
xlim(ax, [-0.7, 12.7]);
ylim(ax, [0, 0.25]);
xticks(ax, 0:3:12);
xlabel(ax, 'i');
ylabel(ax, 'p_i');
legend(ax, 'Location', 'northeast');
savePanel(fig, 'Fig5D');
end

function plotTitrations(ax, titrations, columnPrefix)
%PLOTTITRATIONS  Draw theta(l) for the four enzyme/target conditions of one model.
%   The columns <columnPrefix><condition> of the table hold the curves; they
%   coincide, so line styles as well as colour tell the conditions apart.
conditions = {'1x1x', '2x1x', '1x4x', '2x4x'};
conditionLabels = {'1\times, 1\times', '2\times, 1\times', ...
    '1\times, 4\times', '2\times, 4\times'};
lineStyles = {'-', '--', ':', '-.'};
colors = seriesColors(numel(conditions));
for k = 1:numel(conditions)
    plot(ax, titrations.l, titrations.([columnPrefix, conditions{k}]), lineStyles{k}, ...
        'Color', colors(k, :), 'DisplayName', conditionLabels{k});
end
ax.XScale = 'log';
xlim(ax, [10^-1.2, 10^1.2]);
ylim(ax, [-0.02, 1.02]);
xticks(ax, [0.1, 1, 10]);
xlabel(ax, 'l');
ylabel(ax, '\theta');
legend(ax, 'Location', 'northeast');
end
