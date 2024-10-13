clear all
close all
clc

%% generate packing

% lengths of insertion box in x,y,z
cyl_r = 4e-3;
cyl_h = 1e-2;
cylSize = [cyl_r,cyl_h]; 

% center point location of box
cylCenter = [0,0,cyl_h/2];

% define number of particles and their radii
N_goal = 200;  % goal for number of particles to insert
R_min = 0.37e-3;
R_max = R_min*2;
R_factor = 2*R_min/60e-6;

[sphereCenters, sphereRadii] = placeCylSpheresWeibull(cylSize, cylCenter, N_goal, R_factor, R_min, R_max);

% define properties needed for writing data file 
filename = '../spheres.data';
atom_types = 1;
id = (1:1:N_goal)';
atom_type = ones(N_goal,1);
rho = 1560; % particle density
sphereDensities = rho*ones(N_goal,1);
atoms = horzcat(id,atom_type,2*sphereRadii,sphereDensities,sphereCenters);
simBox = [-4.25e-3 4.25e-3 -4.25e-3 4.25e-3 0 2e-2];
writeFirstDataFile(filename,N_goal,atom_types,simBox,atoms)

%plotWeibullDistribution(sphereRadii,R_factor)

