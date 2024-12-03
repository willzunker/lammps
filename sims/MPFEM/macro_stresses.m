clear all
close all
clc

E = 1e9;
Y = 5e7;
nu = 0.3;
Ro = 0.5;
EY = E/Y;
W = 1.5;

markercolor = {'#fff5cc','#ffe066','#ffcc00','#002db3'};
marker = {'o','s'};
markeredgecolor = '#4d4d4d';
markersize = 13;
linewidth = 0.3;
gcaFontsize = 27;
labelFontsize = 32;
legendFontsize = 27;

%% monodisperse


splice = 200;
data = readmatrix('/Users/willzunker/lammps_mdr_develop/sims/MPFEM/lammps_macro_force_disp.csv');
data = data(1:splice:end, :);
lammps.disp = data(:,1);
lammps.A = (W-2*lammps.disp).^2;
lammps.szz = data(:,2)./lammps.A;
lammps.syy = data(:,3)./lammps.A;
lammps.sxx = data(:,4)./lammps.A;

%splice = 1000;
%data = importdata('fem_macro_force_disp.txt');
%data = data(1:splice:end, :);
%fem.disp = data(:,5);
%fem.A = (W-2*fem.disp).^2;
%fem.szz = data(:,4)./fem.A;
%fem.syy = data(:,3)./fem.A;
%fem.sxx = data(:,2)./fem.A;
%save('fem.mat','fem')
load('fem.mat')

%% differing radii

splice = 200;
data = readmatrix('/Users/willzunker/lammps_mdr_develop/sims/MPFEM/lammps_macro_force_disp_diff_radii.csv');
data = data(1:splice:end, :);
lammps_diff_radii.disp = data(:,1);
lammps_diff_radii.A = (W-2*lammps_diff_radii.disp).^2;
lammps_diff_radii.szz = data(:,2)./lammps_diff_radii.A;
lammps_diff_radii.syy = data(:,3)./lammps_diff_radii.A;
lammps_diff_radii.sxx = data(:,4)./lammps_diff_radii.A;

%splice = 1000;
%data = importdata('fem_macro_force_disp_diff_radii.txt');
%data = data(1:splice:end, :);
%data = data(1:end-63,:);
%fem_diff_radii.disp = data(:,5);
%fem_diff_radii.A = (W-2*fem_diff_radii.disp).^2;
%fem_diff_radii.szz = data(:,4)./fem_diff_radii.A;
%fem_diff_radii.syy = data(:,3)./fem_diff_radii.A;
%fem_diff_radii.sxx = data(:,2)./fem_diff_radii.A;
%save('fem_diff_radii.mat','fem_diff_radii.mat)
load('fem_diff_radii.mat')


%% Combined plots

marker = {'o','s'};
markercolor = {'#819bca','#FFD662'};

figure()
tiledlayout(2,3)

nexttile()
plot(fem.disp./W, fem.sxx./E, marker{1},'MarkerSize', markersize, 'MarkerFaceColor', markercolor{1}, 'MarkerEdgeColor', markeredgecolor, 'LineWidth', linewidth)
hold on
plot(lammps.disp./W, lammps.sxx./E, marker{2},'MarkerSize', markersize, 'MarkerFaceColor', markercolor{2}, 'MarkerEdgeColor', markeredgecolor, 'LineWidth', linewidth)
xlim([0 0.4])
ylim([0 0.65])
set(gcf,'color','w');
set(gca, 'FontSize', gcaFontsize)
set(gca, 'TickLabelInterpreter','latex','XMinorTick','on','YMinorTick','on','Fontsize',gcaFontsize)
xlabel('$\Delta/W$','Interpreter','latex','FontSize', labelFontsize);
ylabel('$\sigma_{xx}/E$','Interpreter','latex','FontSize', labelFontsize);
title('monodisperse','Interpreter','latex','FontSize', legendFontsize)
box on
hl = legend('fem','lammps');
set(hl,'FontSize',legendFontsize,'Interpreter','latex','Location','NorthWest')

nexttile()
plot(fem.disp./W, fem.syy./E, marker{1},'MarkerSize', markersize, 'MarkerFaceColor', markercolor{1}, 'MarkerEdgeColor', markeredgecolor, 'LineWidth', linewidth)
hold on
plot(lammps.disp./W, lammps.syy./E, marker{2},'MarkerSize', markersize, 'MarkerFaceColor', markercolor{2}, 'MarkerEdgeColor', markeredgecolor, 'LineWidth', linewidth)
xlim([0 0.4])
ylim([0 0.65])
set(gcf,'color','w');
set(gca, 'FontSize', gcaFontsize)
set(gca, 'TickLabelInterpreter','latex','XMinorTick','on','YMinorTick','on','Fontsize',gcaFontsize)
xlabel('$\Delta/W$','Interpreter','latex','FontSize', labelFontsize);
ylabel('$\sigma_{yy}/E$','Interpreter','latex','FontSize', labelFontsize);
title('monodisperse','Interpreter','latex','FontSize', legendFontsize)
box on
hl = legend('fem','lammps');
set(hl,'FontSize',legendFontsize,'Interpreter','latex','Location','NorthWest')

nexttile()
plot(fem.disp./W, fem.szz./E, marker{1},'MarkerSize', markersize, 'MarkerFaceColor', markercolor{1}, 'MarkerEdgeColor', markeredgecolor, 'LineWidth', linewidth)
hold on
plot(lammps.disp./W, lammps.szz./E, marker{2},'MarkerSize', markersize, 'MarkerFaceColor', markercolor{2}, 'MarkerEdgeColor', markeredgecolor, 'LineWidth', linewidth)
xlim([0 0.4])
ylim([0 0.65])
set(gcf,'color','w');
set(gca, 'FontSize', gcaFontsize)
set(gca, 'TickLabelInterpreter','latex','XMinorTick','on','YMinorTick','on','Fontsize',gcaFontsize)
xlabel('$\Delta/W$','Interpreter','latex','FontSize', labelFontsize);
ylabel('$\sigma_{zz}/E$','Interpreter','latex','FontSize', labelFontsize);
title('monodisperse','Interpreter','latex','FontSize', legendFontsize)
box on
hl = legend('fem','lammps');
set(hl,'FontSize',legendFontsize,'Interpreter','latex','Location','NorthWest')

nexttile()
plot(fem_diff_radii.disp./W, fem_diff_radii.sxx./E, marker{1},'MarkerSize', markersize, 'MarkerFaceColor', markercolor{1}, 'MarkerEdgeColor', markeredgecolor, 'LineWidth', linewidth)
hold on
plot(lammps_diff_radii.disp./W, lammps_diff_radii.sxx./E, marker{2},'MarkerSize', markersize, 'MarkerFaceColor', markercolor{2}, 'MarkerEdgeColor', markeredgecolor, 'LineWidth', linewidth)
xlim([0 0.4])
ylim([0 0.65])
set(gcf,'color','w');
set(gca, 'FontSize', gcaFontsize)
set(gca, 'TickLabelInterpreter','latex','XMinorTick','on','YMinorTick','on','Fontsize',gcaFontsize)
xlabel('$\Delta/W$','Interpreter','latex','FontSize', labelFontsize);
ylabel('$\sigma_{xx}/E$','Interpreter','latex','FontSize', labelFontsize);
title('tridisperse','Interpreter','latex','FontSize', legendFontsize)
box on
hl = legend('fem','lammps');
set(hl,'FontSize',legendFontsize,'Interpreter','latex','Location','NorthWest')

nexttile()
plot(fem_diff_radii.disp./W, fem_diff_radii.syy./E, marker{1},'MarkerSize', markersize, 'MarkerFaceColor', markercolor{1}, 'MarkerEdgeColor', markeredgecolor, 'LineWidth', linewidth)
hold on
plot(lammps_diff_radii.disp./W, lammps_diff_radii.syy./E, marker{2},'MarkerSize', markersize, 'MarkerFaceColor', markercolor{2}, 'MarkerEdgeColor', markeredgecolor, 'LineWidth', linewidth)
xlim([0 0.4])
ylim([0 0.65])
set(gcf,'color','w');
set(gca, 'FontSize', gcaFontsize)
set(gca, 'TickLabelInterpreter','latex','XMinorTick','on','YMinorTick','on','Fontsize',gcaFontsize)
xlabel('$\Delta/W$','Interpreter','latex','FontSize', labelFontsize);
ylabel('$\sigma_{yy}/E$','Interpreter','latex','FontSize', labelFontsize);
title('tridisperse','Interpreter','latex','FontSize', legendFontsize)
box on
hl = legend('fem','lammps');
set(hl,'FontSize',legendFontsize,'Interpreter','latex','Location','NorthWest')

nexttile()
plot(fem_diff_radii.disp./W, fem_diff_radii.szz./E, marker{1},'MarkerSize', markersize, 'MarkerFaceColor', markercolor{1}, 'MarkerEdgeColor', markeredgecolor, 'LineWidth', linewidth)
hold on
plot(lammps_diff_radii.disp./W, lammps_diff_radii.szz./E, marker{2},'MarkerSize', markersize, 'MarkerFaceColor', markercolor{2}, 'MarkerEdgeColor', markeredgecolor, 'LineWidth', linewidth)
xlim([0 0.4])
ylim([0 0.65])
set(gcf,'color','w');
set(gca, 'FontSize', gcaFontsize)
set(gca, 'TickLabelInterpreter','latex','XMinorTick','on','YMinorTick','on','Fontsize',gcaFontsize)
xlabel('$\Delta/W$','Interpreter','latex','FontSize', labelFontsize);
ylabel('$\sigma_{zz}/E$','Interpreter','latex','FontSize', labelFontsize);
title('tridisperse','Interpreter','latex','FontSize', legendFontsize)
box on
hl = legend('fem','lammps');
set(hl,'FontSize',legendFontsize,'Interpreter','latex','Location','NorthWest')