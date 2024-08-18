clear all
close all
clc

punchRadius = 4; % [mm]
punchArea = pi*punchRadius^2; % [mm^2]

%% Read in lammps data

punchOffset = 0.01;        % [m]
powderFillHeight =   0.0066; %[m]  0.0066; %

lammpsData = readmatrix('upperPunchDispForce.csv');
lammpsData = lammpsData(2:end,:);
lammps.punchPosition = lammpsData(:,1)+punchOffset; % [m]
lammps.strain = (powderFillHeight - lammps.punchPosition)./powderFillHeight; % []
lammps.force = -lammpsData(:,2); % [N]
lammps.punchStress = lammps.force./punchArea; % [MPa]

lammpsData2 = readmatrix('avgStresses.csv');
lammpsData = lammpsData2(1:length(lammpsData),:);
lammps.volFrac = lammpsData(:,4)./(lammps.punchPosition*punchArea*1e-6); % []
lammps.sxx = -lammpsData(:,1)*1e-6.*lammps.volFrac; % [MPa]
lammps.syy = -lammpsData(:,2)*1e-6.*lammps.volFrac; % [MPa]
lammps.szz = -lammpsData(:,3)*1e-6.*lammps.volFrac; % [MPa]
lammps.srr = (lammps.sxx + lammps.syy)./2;          % [MPa]
lammps.szzOsrr = lammps.punchStress./lammps.srr;    % []
lammps.adh_length = lammpsData(:,5);

%% Comparison of experiment to LIGGGHTS
marker = {'o','d','s'};
markercolor = {'#ffe066','#80aaff','#002db3'};
markeredgecolor = '#4d4d4d';
markersize = 10.5;
linewidth = 0.3;
gcaFontsize = 27;
labelFontsize = 32;
legendFontsize = 22;

figure()
tiledlayout(2,2)

% Upper and lower punch stress
nexttile(2)
plot(lammps.strain, lammps.punchStress, marker{3},'MarkerSize', markersize, 'MarkerFaceColor', markercolor{3} , 'MarkerEdgeColor', markeredgecolor, 'LineWidth', linewidth)
ylim([0.001 1.3*max(lammps.punchStress)])
xlim([0 0.5])
set(gcf,'color','w');
set(gca, 'FontSize', gcaFontsize)
set(gca, 'TickLabelInterpreter','latex','XMinorTick','on','YMinorTick','on','Fontsize',gcaFontsize)
xlabel('axial strain []','Interpreter','latex','FontSize', labelFontsize);
ylabel('axial stress [MPa]','Interpreter','latex','FontSize', labelFontsize);
box on
%hl = legend('experiment (upper punch)', 'LIGGGHTS', 'lammps');
hl = legend('lammps');
set(hl,'FontSize',legendFontsize,'Interpreter','latex','Location','NorthWest')

% Radial stress
nexttile(3)
plot(lammps.strain, lammps.srr, marker{3},'MarkerSize', markersize, 'MarkerFaceColor', markercolor{3} , 'MarkerEdgeColor', markeredgecolor, 'LineWidth', linewidth)
ylim([0.001 1.15*max(lammps.srr)])
xlim([0 0.5])
set(gcf,'color','w');
set(gca, 'FontSize', gcaFontsize)
set(gca, 'TickLabelInterpreter','latex','XMinorTick','on','YMinorTick','on','Fontsize',gcaFontsize)
xlabel('axial strain []','Interpreter','latex','FontSize', labelFontsize);
ylabel('radial stress [MPa]','Interpreter','latex','FontSize', labelFontsize);
box on
hl = legend('lammps');
set(hl,'FontSize',legendFontsize,'Interpreter','latex','Location','NorthWest')

% Axial to radial stress ratio
nexttile(4)
plot(lammps.strain, lammps.szzOsrr, marker{3},'MarkerSize', markersize, 'MarkerFaceColor', markercolor{3}, 'MarkerEdgeColor', markeredgecolor, 'LineWidth', linewidth)
ylim([0 5])
xlim([0 0.5])
set(gcf,'color','w');
set(gca, 'FontSize', gcaFontsize)
set(gca, 'TickLabelInterpreter','latex','XMinorTick','on','YMinorTick','on','Fontsize',gcaFontsize)
xlabel('axial strain []','Interpreter','latex','FontSize', labelFontsize);
ylabel('$\sigma_{zz}/\sigma_{rr}$','Interpreter','latex','FontSize', labelFontsize);
box on
hl = legend('lammps');
set(hl,'FontSize',legendFontsize,'Interpreter','latex','Location','NorthEast')

