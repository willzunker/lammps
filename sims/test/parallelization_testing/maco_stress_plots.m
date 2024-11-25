clear all
close all
clc

%% geometric properties of die/fill
punchRadius = 4; % [mm]
punchArea = pi*punchRadius^2; % [mm^2]
punchOffset = 0.01;        % [m]
powderFillHeight =   0.00565; %[m]  0.0066; %

splice = 1;

%% read in serial data
data = readmatrix('/Users/willzunker/lammps/sims/parallelization_testing/E_Y_1_force_disp.csv');
data = data(2:splice:end,:);
serial.punchPosition = data(:,1)+punchOffset; % [m]
serial.strain = (powderFillHeight - serial.punchPosition)./powderFillHeight; % []
serial.force = -data(:,2); % [N]
serial.punchStress = serial.force./punchArea; % [MPa]

data = readmatrix('/Users/willzunker/lammps/sims/parallelization_testing/E_Y_1_avg_stresses.csv');
data = data(2:splice:end,:);
serial.volFrac = data(:,4)./(serial.punchPosition*punchArea*1e-6); % []
serial.sxx = -data(:,1)*1e-6.*serial.volFrac; % [MPa]
serial.syy = -data(:,2)*1e-6.*serial.volFrac; % [MPa]
serial.szz = -data(:,3)*1e-6.*serial.volFrac; % [MPa]
serial.srr = (serial.sxx + serial.syy)./2;          % [MPa]
serial.szzOsrr = serial.punchStress./serial.srr;    % []

%% read in parallel data 1 processors
data = readmatrix('/Users/willzunker/lammps/sims/parallelization_testing/E_Y_1_para1_force_disp.csv');
data = data(2:splice:end,:);
para1.punchPosition = data(:,1)+punchOffset; % [m]
para1.strain = (powderFillHeight - para1.punchPosition)./powderFillHeight; % []
para1.force = -data(:,2); % [N]
para1.punchStress = para1.force./punchArea; % [MPa]

data = readmatrix('/Users/willzunker/lammps/sims/parallelization_testing/E_Y_1_para1_avg_stresses.csv');
data = data(2:splice:end,:);
para1.volFrac = data(:,4)./(para1.punchPosition*punchArea*1e-6); % []
para1.sxx = -data(:,1)*1e-6.*para1.volFrac; % [MPa]
para1.syy = -data(:,2)*1e-6.*para1.volFrac; % [MPa]
para1.szz = -data(:,3)*1e-6.*para1.volFrac; % [MPa]
para1.srr = (para1.sxx + para1.syy)./2;          % [MPa]
para1.szzOsrr = para1.punchStress./para1.srr;    % []

%% read in parallel data 2 processors
data = readmatrix('/Users/willzunker/lammps/sims/parallelization_testing/E_Y_1_para2_force_disp.csv');
data = data(2:splice:end,:);
para2.punchPosition = data(:,1)+punchOffset; % [m]
para2.strain = (powderFillHeight - para2.punchPosition)./powderFillHeight; % []
para2.force = -data(:,2); % [N]
para2.punchStress = para2.force./punchArea; % [MPa]

data = readmatrix('/Users/willzunker/lammps/sims/parallelization_testing/E_Y_1_para2_avg_stresses.csv');
data = data(2:splice:end,:);
para2.volFrac = data(:,4)./(para2.punchPosition*punchArea*1e-6); % []
para2.sxx = -data(:,1)*1e-6.*para2.volFrac; % [MPa]
para2.syy = -data(:,2)*1e-6.*para2.volFrac; % [MPa]
para2.szz = -data(:,3)*1e-6.*para2.volFrac; % [MPa]
para2.srr = (para2.sxx + para2.syy)./2;          % [MPa]
para2.szzOsrr = para2.punchStress./para2.srr;    % []

%% read in parallel data 4 processors
data = readmatrix('/Users/willzunker/lammps/sims/parallelization_testing/E_Y_1_para4_force_disp.csv');
data = data(2:splice:end,:);
para4.punchPosition = data(:,1)+punchOffset; % [m]
para4.strain = (powderFillHeight - para4.punchPosition)./powderFillHeight; % []
para4.force = -data(:,2); % [N]
para4.punchStress = para4.force./punchArea; % [MPa]

data = readmatrix('/Users/willzunker/lammps/sims/parallelization_testing/E_Y_1_para4_avg_stresses.csv');
data = data(2:splice:end,:);
para4.volFrac = data(:,4)./(para4.punchPosition*punchArea*1e-6); % []
para4.sxx = -data(:,1)*1e-6.*para4.volFrac; % [MPa]
para4.syy = -data(:,2)*1e-6.*para4.volFrac; % [MPa]
para4.szz = -data(:,3)*1e-6.*para4.volFrac; % [MPa]
para4.srr = (para4.sxx + para4.syy)./2;          % [MPa]
para4.szzOsrr = para4.punchStress./para4.srr;    % []

%% read in parallel data 8 processors
data = readmatrix('/Users/willzunker/lammps/sims/parallelization_testing/E_Y_1_para8_force_disp.csv');
data = data(2:splice:end,:);
para8.punchPosition = data(:,1)+punchOffset; % [m]
para8.strain = (powderFillHeight - para8.punchPosition)./powderFillHeight; % []
para8.force = -data(:,2); % [N]
para8.punchStress = para8.force./punchArea; % [MPa]

data = readmatrix('/Users/willzunker/lammps/sims/parallelization_testing/E_Y_1_para8_avg_stresses.csv');
data = data(2:splice:end,:);
para8.volFrac = data(:,4)./(para8.punchPosition*punchArea*1e-6); % []
para8.sxx = -data(:,1)*1e-6.*para8.volFrac; % [MPa]
para8.syy = -data(:,2)*1e-6.*para8.volFrac; % [MPa]
para8.szz = -data(:,3)*1e-6.*para8.volFrac; % [MPa]
para8.srr = (para8.sxx + para8.syy)./2;          % [MPa]
para8.szzOsrr = para8.punchStress./para8.srr;    % []

%% plot serial and parallel data together

marker = {'-','-','-','-','-'};
markercolor = {'#ffe066','#80aaff','#1a53ff','#002db3','#002db3'};
markeredgecolor = '#4d4d4d';
markersize = 10.5;
linewidth = 2;
gcaFontsize = 27;
labelFontsize = 32;
legendFontsize = 22;

figure()
tiledlayout(2,2)

% Upper and lower punch stress
nexttile(2)
plot(serial.strain, serial.punchStress, marker{1},'MarkerSize', markersize, 'MarkerFaceColor', markercolor{1} , 'MarkerEdgeColor', markeredgecolor, 'LineWidth', linewidth)
hold on
plot(para1.strain, para1.punchStress, marker{2},'MarkerSize', markersize, 'MarkerFaceColor', markercolor{2} , 'MarkerEdgeColor', markeredgecolor, 'LineWidth', linewidth)
plot(para2.strain, para2.punchStress, marker{3},'MarkerSize', markersize, 'MarkerFaceColor', markercolor{3} , 'MarkerEdgeColor', markeredgecolor, 'LineWidth', linewidth)
plot(para4.strain, para4.punchStress, marker{4},'MarkerSize', markersize, 'MarkerFaceColor', markercolor{4} , 'MarkerEdgeColor', markeredgecolor, 'LineWidth', linewidth)
plot(para8.strain, para8.punchStress, marker{5},'MarkerSize', markersize, 'MarkerFaceColor', markercolor{5} , 'MarkerEdgeColor', markeredgecolor, 'LineWidth', linewidth)
ylim([0 1.3*max(serial.punchStress)])
xlim([0 0.5])
set(gcf,'color','w');
set(gca, 'FontSize', gcaFontsize)
set(gca, 'TickLabelInterpreter','latex','XMinorTick','on','YMinorTick','on','Fontsize',gcaFontsize)
xlabel('axial strain []','Interpreter','latex','FontSize', labelFontsize);
ylabel('axial stress [MPa]','Interpreter','latex','FontSize', labelFontsize);
box on
%hl = legend('experiment (upper punch)', 'LIGGGHTS', 'lammps');
hl = legend('serial', 'parallel: 1 proc.', 'parallel: 2 proc.', 'parallel: 4 proc.','parallel: 8 proc.');
set(hl,'FontSize',legendFontsize,'Interpreter','latex','Location','NorthWest')

% Radial stress
nexttile(3)
plot(serial.strain, serial.srr, marker{1},'MarkerSize', markersize, 'MarkerFaceColor', markercolor{1} , 'MarkerEdgeColor', markeredgecolor, 'LineWidth', linewidth)
hold on
plot(para1.strain, para1.srr, marker{2},'MarkerSize', markersize, 'MarkerFaceColor', markercolor{2} , 'MarkerEdgeColor', markeredgecolor, 'LineWidth', linewidth)
plot(para2.strain, para2.srr, marker{3},'MarkerSize', markersize, 'MarkerFaceColor', markercolor{3} , 'MarkerEdgeColor', markeredgecolor, 'LineWidth', linewidth)
plot(para4.strain, para4.srr, marker{4},'MarkerSize', markersize, 'MarkerFaceColor', markercolor{4} , 'MarkerEdgeColor', markeredgecolor, 'LineWidth', linewidth)
plot(para8.strain, para8.srr, marker{5},'MarkerSize', markersize, 'MarkerFaceColor', markercolor{5} , 'MarkerEdgeColor', markeredgecolor, 'LineWidth', linewidth)
ylim([0 1.15*max(serial.srr)])
xlim([0 0.5])
set(gcf,'color','w');
set(gca, 'FontSize', gcaFontsize)
set(gca, 'TickLabelInterpreter','latex','XMinorTick','on','YMinorTick','on','Fontsize',gcaFontsize)
xlabel('axial strain []','Interpreter','latex','FontSize', labelFontsize);
ylabel('radial stress [MPa]','Interpreter','latex','FontSize', labelFontsize);
box on
%hl = legend('experimental (sensor average)', 'LIGGGHTS', 'lammps');
hl = legend('serial', 'parallel: 1 proc.', 'parallel: 2 proc.', 'parallel: 4 proc.','parallel: 8 proc.');
set(hl,'FontSize',legendFontsize,'Interpreter','latex','Location','NorthWest')

% Axial to radial stress ratio
nexttile(4)
plot(serial.strain, serial.szzOsrr, marker{1},'MarkerSize', markersize, 'MarkerFaceColor', markercolor{1}, 'MarkerEdgeColor', markeredgecolor, 'LineWidth', linewidth)
hold on
plot(para1.strain, para1.szzOsrr, marker{2},'MarkerSize', markersize, 'MarkerFaceColor', markercolor{2}, 'MarkerEdgeColor', markeredgecolor, 'LineWidth', linewidth)
plot(para2.strain, para2.szzOsrr, marker{3},'MarkerSize', markersize, 'MarkerFaceColor', markercolor{3}, 'MarkerEdgeColor', markeredgecolor, 'LineWidth', linewidth)
plot(para4.strain, para4.szzOsrr, marker{4},'MarkerSize', markersize, 'MarkerFaceColor', markercolor{4}, 'MarkerEdgeColor', markeredgecolor, 'LineWidth', linewidth)
plot(para8.strain, para8.szzOsrr, marker{5},'MarkerSize', markersize, 'MarkerFaceColor', markercolor{5}, 'MarkerEdgeColor', markeredgecolor, 'LineWidth', linewidth)
ylim([0 5])
xlim([0 0.5])
set(gcf,'color','w');
set(gca, 'FontSize', gcaFontsize)
set(gca, 'TickLabelInterpreter','latex','XMinorTick','on','YMinorTick','on','Fontsize',gcaFontsize)
xlabel('axial strain []','Interpreter','latex','FontSize', labelFontsize);
ylabel('$\sigma_{zz}/\sigma_{rr}$','Interpreter','latex','FontSize', labelFontsize);
box on
hl = legend('serial', 'parallel: 1 proc.', 'parallel: 2 proc.', 'parallel: 4 proc.','parallel: 8 proc.');
set(hl,'FontSize',legendFontsize,'Interpreter','latex','Location','NorthEast')
