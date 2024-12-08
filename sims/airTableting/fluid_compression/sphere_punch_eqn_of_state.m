clear all
close all
clc

data = readmatrix('sphere_punch_eqn_of_state.csv');
volume = data(2:end,1)*8e-3*8e-3;
pressure = data(2:end,2);
volume = mean(reshape(volume,200,[]),1);
pressure = mean(reshape(pressure,200,[]),1);
N = data(1,3);
T = data(1,4);
kb = 1.380649e-23;
Vo = 8e-3*8e-3*1e-2;

gcaFontsize = 27;
labelFontsize = 32;
legendFontsize = 22;

pressure_ana = (1./volume)*N*T*kb;


%%
p = 101325;
V = 1;
R = 8.314;
T_air = 293;

n = p*V/(R*T_air);

n_air = n*Vo;

pressure_air = (1./volume)*n_air*R*T_air;

figure()
plot(volume/Vo,pressure,'LineWidth',3)
hold on
plot(volume/Vo,pressure_ana,'LineWidth',3)
plot(volume/Vo,pressure_air,'LineWidth',3)
set(gcf,'color','w');
set(gca, 'FontSize', gcaFontsize)
set(gca, 'TickLabelInterpreter','latex','XMinorTick','on','YMinorTick','on','Fontsize',gcaFontsize)
xlabel('$V/V_o$','Interpreter','latex','FontSize', labelFontsize);
ylabel('$P$ [Pa]','Interpreter','latex','FontSize', labelFontsize);
box on
hl = legend('lammps', '$P = \frac{Nk_bT}{V}$', 'air');
set(hl,'FontSize',legendFontsize,'Interpreter','latex','Location','NorthEast')