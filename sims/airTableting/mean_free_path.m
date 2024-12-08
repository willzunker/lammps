clear all
close all
clc

%% mean free path

N_reset = 20; % check
dt = 1e-6; % check
dt_srd = N_reset*dt;
kb = 1.380649e-23;
T_mpcd = 9.4e15; % check
R_big = 0.4e-3; 
rho_big = 1560;
m_big = 4/3*pi*R_big^3*rho_big;
m_mpcd = m_big*1e-4; % check
lambda = dt_srd*sqrt(kb*T_mpcd/m_mpcd)

N_mpcd = 500000; % check

ideal_constant_lammps = N_mpcd*kb*T_mpcd

%% pV = nRT

p_air = 101325;
V_unit = 1;
R = 8.314;
T_air = 293;

n_unit = p_air*V_unit/(R*T_air);
V_die = 8e-3*8e-3*1e-2;
n_air = n_unit*V_die;

ideal_constant_air = n_air*R*T_air

%% settling velocity

m_air = 1.29*V_die;
m_mpcd_all = m_mpcd*N_mpcd;
rho_mpcd = m_mpcd_all/V_die;

rho_air = 1.29;
R_avicel = 50e-6;
g = 9.81;
nu_air = 1.48e-5;

% Fg = 4*pi*(rho_m-rho_w)*g*R^3/3;
% Fd = 6*pi*nu*rho_w*R*v;

v_avicel = 2/9*rho_big*g*R_avicel^2/(nu_air*rho_air)
v_avicel = 2/9*rho_big*g*R_big^2/(nu_air*rho_air)
g_mpcd = 25*g; % check
v_big = 2/9*rho_big*g_mpcd*R_big^2/(nu_air*rho_mpcd)

%% viscocity
a = 2*R_big/4; % check
M = N_mpcd/(V_die/a^3)
nu_air 
nu_mpcd = a^2/(18*dt_srd)*(1-(1-exp(-M))./M) + kb*T_mpcd*dt_srd*(M+2)./(4*m_mpcd*(M-1))

term1 = a^2/(18*dt_srd)*(1-(1-exp(-M))./M)
term2 = kb*T_mpcd*dt_srd*(M+2)./(4*m_mpcd*(M-1))

M = (0:0.001:5)';
c = 1/18*(1-(1-exp(-M))./M) + 1/4*(M+2)./(M-1)

figure()
plot(M,c)
ylim([-20 20])

nu = a^2./(18*dt_srd).*(1-(1-exp(-M))./M) + kb*T_mpcd*dt_srd*(M+2)./(4*m_mpcd*(M-1)); 

figure()
semilogx(M,nu)
ylim([-1e-3 1e-3])

%%

compression_steps = 10000;
zhi = 10;
dt = 1e-4;
V_punch = (1/10000)*0.25*zhi/1e-3


E = 1e4;
R = 0.5e-3;
k = E*R;
rho = 1560;
V = (4/3)*pi*R^3;
m = rho*V;

dt_crit = sqrt(m/k)


R_dan = 0.5;
rho_dan = 1.980;
V_dan = (4/3)*pi*R_dan^3;
m_dan = rho_dan*V_dan

%% viscocity
a = 2e-4;
dt = 1e-5;
N_reset = 20;
dt_srd = N_reset*dt;
M = (0:0.001:5)';
kb = 1.380649e-23;
T_mpcd = 1e13;
R_big = 0.5e-3;
rho_big = 1560;
mf = 4/3*pi*R_big^3*rho_big*1e-3;
nu_air = a^2/(18*dt_srd)*(1-(1-exp(-M))./M) + kb*T_mpcd*dt_srd*(M+2)./(4*mf*(M-1));

figure()
semilogx(M,nu_air)
