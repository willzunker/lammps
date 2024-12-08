clear all
%close all
clc

% %% Read in Vertex data
% rawdata = readmatrix('CompactionSimulatorData_Avicel102/19574_MV_12.5kN_12.5_-500KNs_2s_-12.75mm_-5mm_LK_25kHz_Avicelph102-P216829735_C01R02T01L01.txt');
% %rawdata = readmatrix('CompactionSimulatorData_Avicel102/19574_MV_12.5kN_12.5_-500KNs_2s_-12.75mm_-5mm_LK_25kHz_Avicelph102-P216829735_C02R02T01L01.txt');
% punchRadius = 4; % [mm]
% punchArea = pi*punchRadius^2; % [mm^2]
% initialFillPos = -6.9; % [mm]
% 
% vertex.time = rawdata(:,1); % [sec]
% 
% vertex.upperPunchPos = rawdata(:,2); % [mm]
% vertex.lowerPunchPos = rawdata(:,3); % [mm]
% vertex.powderStrain = (initialFillPos - vertex.upperPunchPos)./(initialFillPos - mean(vertex.lowerPunchPos(1:50)));
% vertex.upperPunchForce = rawdata(:,4); % [kN]
% vertex.lowerPunchForce = rawdata(:,5); % [kN]
% vertex.upperPunchStress = vertex.upperPunchForce*10^3./punchArea; % [MPa]
% vertex.lowerPunchStress = vertex.lowerPunchForce*10^3./punchArea; % [MPa]
% vertex.avgPunchStress = (vertex.upperPunchStress + vertex.lowerPunchStress)./2; % [MPa]
% vertex.upperRadialStress = rawdata(:,8); % [MPa]
% vertex.lowerRadialStress = rawdata(:,10); % [MPa]
% vertex.avgRadialStress = (vertex.upperRadialStress + vertex.lowerRadialStress)./2; % [MPa]
% vertex.szzOsrr = vertex.upperPunchStress./vertex.avgRadialStress; % []

% %% Read in LIGGGHTS data
% punchOffset = 0.01;        % [m]
% powderFillHeight = 0.0062; % [m]
% 
% liggghtsData = readmatrix('../benchmarkLIGGGHTsMDRmodel/benchmarkMDRmodelCSVs/punchForce.csv');
% liggghts.punchPosition = liggghtsData(:,1)+punchOffset; % [m]
% liggghts.strain = (powderFillHeight - liggghts.punchPosition)./powderFillHeight; % []
% liggghts.force = liggghtsData(:,2); % [N]
% liggghts.punchStress = liggghts.force./punchArea; % [MPa]
% 
% liggghtsData = readmatrix('../benchmarkLIGGGHTsMDRmodel/benchmarkMDRmodelCSVs/atomStresses.csv');
% liggghts.volFrac = liggghtsData(:,4)./(liggghts.punchPosition*punchArea*1e-6); % []
% liggghts.sxx = -liggghtsData(:,1)*1e-6.*liggghts.volFrac; % [MPa]
% liggghts.syy = -liggghtsData(:,2)*1e-6.*liggghts.volFrac; % [MPa]
% liggghts.szz = -liggghtsData(:,3)*1e-6.*liggghts.volFrac; % [MPa]
% liggghts.srr = (liggghts.sxx + liggghts.syy)./2;          % [MPa]
% liggghts.szzOsrr = liggghts.punchStress./liggghts.srr;    % []

%% Read in lammps data

punchOffset = 0.01;        % [m]
powderFillHeight = 0.006; % [m]
punchArea = 0.8e-3^2;

lammpsData = readmatrix('force_disp_no_air.csv');
lammpsData = lammpsData(2:end,:);
lammps.punchPosition = lammpsData(:,1)+punchOffset; % [m]
lammps.strain = (powderFillHeight - lammps.punchPosition)./powderFillHeight; % []
lammps.force = -lammpsData(:,2); % [N]
lammps.punchStress = lammps.force./punchArea; % [MPa]


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
%plot(vertex.powderStrain, vertex.upperPunchStress, marker{1},'MarkerSize', markersize, 'MarkerFaceColor', markercolor{1}, 'MarkerEdgeColor', markeredgecolor, 'LineWidth', linewidth)
hold on
plot(lammps.strain, lammps.punchStress, marker{2},'MarkerSize', markersize, 'MarkerFaceColor', markercolor{2} , 'MarkerEdgeColor', markeredgecolor, 'LineWidth', linewidth)
%ylim([0.001 1.15*max(vertex.avgPunchStress)])
xlim([0 0.5])
set(gcf,'color','w');
set(gca, 'FontSize', gcaFontsize)
set(gca, 'TickLabelInterpreter','latex','XMinorTick','on','YMinorTick','on','Fontsize',gcaFontsize)
xlabel('axial strain []','Interpreter','latex','FontSize', labelFontsize);
ylabel('axial stress [MPa]','Interpreter','latex','FontSize', labelFontsize);
box on
hl = legend('lammps');
set(hl,'FontSize',legendFontsize,'Interpreter','latex','Location','NorthWest')

% % Radial stress
% nexttile(3)
% plot(vertex.powderStrain, vertex.avgRadialStress, marker{1},'MarkerSize', markersize, 'MarkerFaceColor', markercolor{1} , 'MarkerEdgeColor', markeredgecolor, 'LineWidth', linewidth)
% hold on
% plot(liggghts.strain, liggghts.srr, marker{2},'MarkerSize', markersize, 'MarkerFaceColor', markercolor{2} , 'MarkerEdgeColor', markeredgecolor, 'LineWidth', linewidth)
% plot(lammps.strain, lammps.srr, marker{3},'MarkerSize', markersize, 'MarkerFaceColor', markercolor{3} , 'MarkerEdgeColor', markeredgecolor, 'LineWidth', linewidth)
% ylim([0.001 1.15*max(vertex.avgRadialStress)])
% xlim([0 0.5])
% set(gcf,'color','w');
% set(gca, 'FontSize', gcaFontsize)
% set(gca, 'TickLabelInterpreter','latex','XMinorTick','on','YMinorTick','on','Fontsize',gcaFontsize)
% xlabel('axial strain []','Interpreter','latex','FontSize', labelFontsize);
% ylabel('radial stress [MPa]','Interpreter','latex','FontSize', labelFontsize);
% box on
% hl = legend('experimental (sensor average)', 'LIGGGHTS', 'lammps');
% set(hl,'FontSize',legendFontsize,'Interpreter','latex','Location','NorthWest')
% 
% % Axial to radial stress ratio
% nexttile(4)
% hold on
% plot(vertex.powderStrain, vertex.szzOsrr, marker{1},'MarkerSize', markersize, 'MarkerFaceColor', markercolor{1}, 'MarkerEdgeColor', markeredgecolor, 'LineWidth', linewidth)
% plot(liggghts.strain, liggghts.szzOsrr, marker{2},'MarkerSize', markersize, 'MarkerFaceColor', markercolor{2}, 'MarkerEdgeColor', markeredgecolor, 'LineWidth', linewidth)
% plot(lammps.strain, lammps.szzOsrr, marker{3},'MarkerSize', markersize, 'MarkerFaceColor', markercolor{3}, 'MarkerEdgeColor', markeredgecolor, 'LineWidth', linewidth)
% ylim([0 5])
% xlim([0 0.5])
% set(gcf,'color','w');
% set(gca, 'FontSize', gcaFontsize)
% set(gca, 'TickLabelInterpreter','latex','XMinorTick','on','YMinorTick','on','Fontsize',gcaFontsize)
% xlabel('axial strain []','Interpreter','latex','FontSize', labelFontsize);
% ylabel('$\sigma_{zz}/\sigma_{rr}$','Interpreter','latex','FontSize', labelFontsize);
% box on
% hl = legend('experiment', 'LIGGGHTS', 'lammps');
% set(hl,'FontSize',legendFontsize,'Interpreter','latex','Location','NorthEast')
% 
% 
% % Set the figure size to match the paper size
% fig = gcf; % Get current figure handle
% fig.Units = 'inches'; % Set units to inches
% fig.PaperSize = [24 16];
% fig.Position = [0, 0, 24, 16]; % Set figure position and size to match paper size
% 
% % Save the figure as a PDF
% print('/Users/willzunker/GradSchool/Conferences/Pfizer_talk/LIGGGHTS_experiment_comparison.pdf', '-dpdf', '-fillpage');



















%% plot axial strain vs steps

% figure()
% steps = (1:1:150000)';
% 
% plot(steps,liggghtsAdhOff.strain, '<','MarkerSize', markersize, 'MarkerFaceColor', markercolor{1}, 'MarkerEdgeColor', markeredgecolor, 'LineWidth', linewidth)
% hold on
% %ylim([0 0.4])
% %xlim([0 300])
% set(gcf,'color','w');
% set(gca, 'FontSize', gcaFontsize)
% set(gca, 'TickLabelInterpreter','latex','XMinorTick','on','YMinorTick','on','Fontsize',gcaFontsize)
% xlabel('steps','Interpreter','latex','FontSize', labelFontsize);
% ylabel('axial strain []','Interpreter','latex','FontSize', labelFontsize);
% box on
% set(hl,'FontSize',legendFontsize,'Interpreter','latex','Location','NorthWest')
