clear all
close all
clc

%% pV = nRT

p_air = 101325;
V_unit = 1;
R = 8.314;
T_air = 293;

n_unit = p_air*V_unit/(R*T_air);
V_die = 8e-3*8e-3*1e-2;
n_air = n_unit*V_die;

ideal_constant_air = n_air*R*T_air;

%% mean free path

N_mpcd = 3e5; % check
N_reset = 25; % check
dt = 2e-07; % check
dt_srd = N_reset*dt;
kb = 1.380649e-23;
T_mpcd = ideal_constant_air/(N_mpcd*kb); % check
R_big = 0.4e-3; 
rho_big = 1560;
rho_big_heavy = 100*rho_big;
m_big = 4/3*pi*R_big^3*rho_big;
m_mpcd = m_big*3e-4; % check
lambda = dt_srd*sqrt(kb*T_mpcd/m_mpcd);
ideal_constant_mpcd = N_mpcd*kb*T_mpcd;

fprintf('lambda is: %s\n', lambda);
fprintf('ideal gas constant (air, mpcd): %f, %f. T_mpcd = %s\n', ideal_constant_air, ideal_constant_mpcd,T_mpcd);

%% viscocity
nu_air = 1.48e-5;

a = 2*R_big/4; % check
M = N_mpcd/(V_die/a^3);
nu_mpcd = a^2/(18*dt_srd)*(1-(1-exp(-M))./M) + lambda^2*(M+2)./(4*dt_srd*(M-1));

fprintf('viscosity (air, mpcd): %s, %s\n', nu_air, nu_mpcd);

fprintf('M: %f \n', M);

fprintf('lambda/a: %f\n', lambda/a);

%% settling velocity

m_air = 1.29*V_die;
m_mpcd_all = m_mpcd*N_mpcd;
rho_mpcd = m_mpcd_all/V_die;

fprintf('Big to mpcd density ratio: %f. mass factor %s \n', rho_big/rho_mpcd, m_mpcd/m_big);
fprintf('Big heavy to mpcd density ratio: %f. \n', rho_big_heavy/rho_mpcd);

rho_air = 1.29;
R_avicel = 50e-6;
g = 9.81;

% Fg = 4*pi*(rho_m-rho_w)*g*R^3/3;
% Fd = 6*pi*nu*rho_w*R*v;

v_avicel_air = 2/9*rho_big*g*R_avicel^2/(nu_air*rho_air);
v_big_air = 2/9*rho_big*g*R_big^2/(nu_air*rho_air);
v_big_mpcd = 2/9*rho_big*g*R_big^2/(nu_mpcd*rho_mpcd);
v_big_heavy_mpcd = 2/9*rho_big_heavy*g*R_big^2/(nu_mpcd*rho_mpcd);
v_heavy_target = 2;
g_factor_avicel = v_avicel_air/v_big_mpcd; % check
g_factor_big = v_big_air/v_big_mpcd;
g_factor_big_heavy = v_heavy_target/v_big_heavy_mpcd;

fprintf('settling velocity avicel sized particle (air, mpcd): %f, %f. gravity factor required %f\n', v_avicel_air, v_big_mpcd, g_factor_avicel);
fprintf('settling velocity big sized particle (air, mpcd): %f, %f. gravity factor required %f\n', v_big_air, v_big_mpcd, g_factor_big);
fprintf('settling velocity big heavy sized particle (target, mpcd): %f, %f. gravity factor required %f\n', v_heavy_target, v_big_heavy_mpcd, g_factor_big_heavy);

%% volumetric flow rate

R_die = 4e-3;
w_die = 50e-6;
v_punch = 0.25;
Q = v_punch*2*pi*R_die*w_die;

dx_avg = 0.0031;

wOnu = w_die^3/(nu_air);
wOnudx = w_die^3/(nu_air*dx_avg);
w_required = (wOnu*nu_mpcd)^(1/3);
w_required_tight = (wOnudx*nu_mpcd*R_big)^(1/3);

fprintf('die-punch gap : %s, %s, %s\n', w_die, w_required,w_required_tight);

%Re_piston = v_punch*w/nu_air;
%Re_particle = v_punch*R_big/nu_air;
%v_punch_required = Re_particle*nu_mpcd/R_big;
%w_required = Re_piston*nu_mpcd/v_punch_required;
%fprintf('To match the Re the piston speed and gap must be : %f, %f\n', v_punch_required, w_required);

%% Carman-Kozeny equation

porosity_avicel = 0.8;
kappa_avicel = 1^2*porosity_avicel^3*R_avicel^2/(180*(1-porosity_avicel)^2);
kOmu_air = kappa_avicel/(rho_air*nu_air);

powder_fill_height = 0.0066;
V_total = pi*R_die^2*powder_fill_height; 
V_big_total = 4/3*pi*R_big^3*700;
porosity_big = 1 - V_big_total/V_total;
kappa_big =  1^2*porosity_big^3*R_big^2/(180*(1-porosity_big)^2);
kOmu_mpcd = kappa_big/(rho_mpcd*nu_mpcd);
fprintf('permiability to viscosity ratio (air, mpcd): %s, %s\n', kOmu_air,kOmu_mpcd);

%% error checks
if M < 1
    error("M < 1, increase a or N_mpcd")
end

if lambda/a < 0.6 
    error("lambda/a < 0.6, shifts will be required")
end

if lambda > R_big
    error("lambda > R_big, collisions between mpcd and big particles might not occur")
end

if rho_mpcd > rho_big
    error("rho_mpcd > rho_big, density of mpcd is larger than big particles issues may occur")
end

%% plot settling speed 

velocity_time = readmatrix('fluid_compression/particle_settle_velocity_time.csv');

t_lammps = velocity_time(:,1);
v_lammps = velocity_time(:,2);

n = length(velocity_time);
new_length = floor(n / 400) * 400;

t_lammps = t_lammps(1:new_length);
v_lammps = v_lammps(1:new_length);

t_lammps = mean(reshape(t_lammps,400,[]),1);
v_lammps = mean(reshape(v_lammps,400,[]),1);

t_ana = t_lammps';

c = 9/2*nu_mpcd*rho_mpcd/(R_big^2*rho_big);
g = v_big_air/v_big_mpcd*9.81;
v_ana = g/c*(exp(-c*t_ana) - 1);

gcaFontsize = 27;
labelFontsize = 32;
legendFontsize = 22;

figure()
plot(t_lammps,v_lammps,'-sk','MarkerSize',8,'LineWidth',2)
hold on
plot(t_ana,v_ana,'-b','LineWidth',2)
set(gcf,'color','w');
set(gca, 'FontSize', gcaFontsize)
set(gca, 'TickLabelInterpreter','latex','XMinorTick','on','YMinorTick','on','Fontsize',gcaFontsize)
xlabel('$t$ [s]','Interpreter','latex','FontSize', labelFontsize);
ylabel('$v$ [m/s]','Interpreter','latex','FontSize', labelFontsize);
box on
hl = legend('lammps', 'analytical');
set(hl,'FontSize',legendFontsize,'Interpreter','latex','Location','NorthEast')


%% print off 
%syms dt_ m_
%nu_air = 1e-3;
%dt_required = vpasolve( a^2/(18*dt_)*(1-(1-exp(-M))/M) + kb*T_mpcd*dt_*(M+2)/(4*m_mpcd*(M-1)) == nu_air, dt_, 1 )
%m_required = vpasolve( a^2/(18*dt_srd)*(1-(1-exp(-M))/M) + kb*T_mpcd*dt_srd*(M+2)/(4*m_*(M-1)) == nu_air, m_, 1 )
