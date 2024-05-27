clear all
close all
clc

data = readmatrix('eqn_of_state.csv');
volume = data(:,1)*10*10;
pressure = data(:,2);
N = data(1,3);
T = data(1,4);
kb = 1.380649e-23;
Vo = 10*10*10;

Ao = 10*10;
A = data(:,1)*10;
A_adjust = (6*Ao)./(Ao*2 + 4*A);

pressure = pressure.*A_adjust;

gcaFontsize = 27;
labelFontsize = 32;
legendFontsize = 22;

pressure_ana = (1./volume)*N*T*kb;

figure()
plot(volume/Vo,pressure,'LineWidth',3)
hold on
plot(volume/Vo,pressure_ana,'LineWidth',3)
set(gcf,'color','w');
set(gca, 'FontSize', gcaFontsize)
set(gca, 'TickLabelInterpreter','latex','XMinorTick','on','YMinorTick','on','Fontsize',gcaFontsize)
xlabel('$V/V_o$','Interpreter','latex','FontSize', labelFontsize);
ylabel('$P$ [Pa]','Interpreter','latex','FontSize', labelFontsize);
box on
hl = legend('lammps', '$P = \frac{Nk_bT}{V}$');
set(hl,'FontSize',legendFontsize,'Interpreter','latex','Location','NorthEast')