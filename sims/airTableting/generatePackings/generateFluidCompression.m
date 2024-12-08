clear all
close all
clc

%% generate square punch

% lengths of insertion box in x,y,z
box_x = 8e-3;
box_y = 8e-3;
box_z = 0.99e-2;
boxSize = [box_x, box_y, box_z]; 

% center point location of box
boxCenter = [0,0,box_z/2];

% define number of particles and their radii
N_goal = 0;  % goal for number of particles to insert

R_punch = box_x/20;
N_x = 30;
N_y = 30; 
atom_types = 2;

[sphereCenters, sphereRadii] = generateSquarePunch(R_punch,boxSize,boxCenter,N_x,N_y);

% define properties needed for writing data file 
filename = '../fluid_compression/sphere_punch.data';
N_punchMax = length(sphereRadii);
id = ((N_goal+1):1:(N_goal+N_punchMax))';
N_punch = length(id);
atom_type = 1*ones(N_punch,1);
rho = 1560; % particle density
sphereDensities = rho*ones(N_punch,1);
simBox = [-4e-3 4e-3 -4e-3 4e-3 0 1e-2];
atoms = horzcat(id,atom_type,2*sphereRadii,sphereDensities,sphereCenters);
writeFirstDataFile(filename,N_punch,atom_types,simBox,atoms)



