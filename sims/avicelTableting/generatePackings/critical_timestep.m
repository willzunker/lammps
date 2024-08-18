clear all
close all
clc

% timestep

rho = 1560;
R = 0.4e-3;
E = 7.5e9;

dt = sqrt(rho*R^2/E)

E = 1e10;
nu = 0.3

Reff = (1/0.65e-3 + 1/0.65e-3)^(-1)
Ecom = E/(1-nu^2)

kt = Reff*Ecom*2/7


% volume
xl = -4e-3;
xh = 4e-3;
yl = -4e-3;
yh = 4e-3;
zl = 0;
zh = 1.2e-2;
V1 = (xh-xl)*(yh-yl)*(zh-zl)

zl_f = 1e-2;
zh_f = 2e-2;
c = sqrt(V1/((xh-xl)*(yh-yl)*(zh_f-zl_f)))
xl_f = c*xl;
xh_f = c*xh;
yl_f = c*yl;
yh_f = c*xh;

V2 = (xh_f-xl_f)*(yh_f-yl_f)*(zh_f-zl_f)

