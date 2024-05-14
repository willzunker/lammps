clear all
close all
clc

N = 200;
dt = 1e-5;
kb = 1.380649e-23;
T = 1e11;
R_big = 0.5e-3;
rho_big = 1560;
m = 4/3*pi*R_big^3*rho_big*1e-3;
lambda = N*dt*sqrt(kb*T/m)




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


