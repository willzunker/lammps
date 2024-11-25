function plot_residual_stresses_strain(results_file,Y,strain)

for i = 1 : length(strain)
    force_disp_file = strcat('/Users/willzunker/lammps/sims/residual_stress_study/',results_file,'/',results_file,'_',num2str(i),'_force_disp.csv');
    avg_stresses_file = strcat('/Users/willzunker/lammps/sims/residual_stress_study/',results_file,'/',results_file,'_',num2str(i),'_avg_stresses.csv');
    force_disp = abs(readmatrix(force_disp_file));
    avg_stresses = abs(readmatrix(avg_stresses_file));
    [max_force,max_index] = max(force_disp(:,2));
    threshold = 0.01*max_force;
    indices = find(force_disp(max_index:end,2) < threshold);
    residual_index = indices(1)+(max_index-1);
    residual_stress(i) = mean(avg_stresses(residual_index,1:2));
end

marker = {'o'};
markercolor = {'#ffe066'};
markeredgecolor = '#4d4d4d';
markersize = 10.5;
linewidth = 0.3;
gcaFontsize = 27;
labelFontsize = 32;
legendFontsize = 22;

figure()
plot(strain', residual_stress./Y, marker{1},'MarkerSize', markersize, 'MarkerFaceColor', markercolor{1} , 'MarkerEdgeColor', markeredgecolor, 'LineWidth', linewidth)
ylim([0 1.15*max(residual_stress./Y)])
%xlim([0 0.5])
set(gcf,'color','w');
set(gca, 'FontSize', gcaFontsize)
set(gca, 'TickLabelInterpreter','latex','XMinorTick','on','YMinorTick','on','Fontsize',gcaFontsize)
xlabel('$\epsilon$','Interpreter','latex','FontSize', labelFontsize);
ylabel('$\sigma_\textrm{res}/Y$','Interpreter','latex','FontSize', labelFontsize);
box on
end