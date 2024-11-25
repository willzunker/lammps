clear all
close all
clc

% geometric and material inputs
rho = 1560;   % density
R = 0.0961e-3;  % particle radius
E = 5e9;      % Young's modulus
nu = 0.4;     % Poisson's ratio
Y = 1.5e6;    % Yield stress

% calculated values from inputs 
V = 4/3*pi*R^3;         % particle volume
m = rho*V;              % particle mass
Eeff = E/(1-nu^2);      % composite plane strain modulus
kappa = E/(3*(1-2*nu)); % bulk modulus

% assume worst case overlap to find effective radius
deltaOR = 0.75;
py = Y*(1.75*exp(-4.4*deltaOR) + 1);
amax = deltaOR*R;
A = 4*py/Eeff*amax;
B = 2*amax;
Reff = (B/2)^2/(A/2);

% different stiffness measures
k1 = E*R; % scaling based on initial values
k2 = 1/V*(4*pi*R^2)*kappa*R^2; % bulk stiffness
k3 = 2*Eeff*Reff; % mdr unloading stiffness based on updated effective radius

% critical timestep
dt1 = sqrt(m/k1)
dt2 = sqrt(m/k2)
dt3 = sqrt(m/k3)

%dt = sqrt(rho*R^2/(E))
%dt = 0.8e-7;
%punchFactor = 0.6915;
%max_compression_step = round(0.02/dt) + round(punchFactor*1e-2/0.25/dt)
%neighbor = 1.75*R*1.5

% relaxation time of damper

CoR = 0.5;
beta = -log(CoR)/(sqrt(log(CoR)^2 + pi^2));
Eeff = E/(1-nu^2);
tau = beta*R*sqrt(rho*Eeff)/Eeff




