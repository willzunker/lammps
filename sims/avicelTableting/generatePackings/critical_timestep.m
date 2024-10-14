clear all
close all
clc

% timestep

rho = 1560;
R = 0.133e-3;
E = 9e9;
V = 4/3*pi*R^3;
m = rho*V;
nu = 0.3;
kappa = E/(3*(1-2*nu));
k = 1/V*(4*pi*R^2)*kappa*R^2;

dt = sqrt(rho*R^2/E)

dt = sqrt(m/k)

dt = 0.8e-7;
punchFactor = 0.6915;
max_compression_step = round(0.02/dt) + round(punchFactor*1e-2/0.25/dt)

neighbor = 1.75*R*1.5

