clear all
close all
clc

%% generate packing

% lengths of insertion box in x,y,z
box_x = 8e-3;
box_y = 8e-3;
box_z = 1e-2*0.975;
boxSize = [box_x, box_y, box_z]; 

% center point location of box
boxCenter = [0,0,box_z/2];

% define number of particles and their radii
N_goal = 3000;  % goal for number of particles to insert
R_mean = 0.24e-3;%0.4e-3;   % average sphere radius
R_std = R_mean/20;   % standard deviation of sphere radius

[sphereCenters, sphereRadii] = placeBoxSpheres(boxSize, boxCenter, N_goal, R_mean, R_std);

% define properties needed for writing data file 
filename = '../spheres.data';
atom_types = 3;
id = (1:1:N_goal)';
atom_type = ones(N_goal,1);
rho = 1560; % particle density
sphereDensities = rho*ones(N_goal,1);
atoms = horzcat(id,atom_type,2*sphereRadii,sphereDensities,sphereCenters);
simBox = [-5e-3 5e-3 -5e-3 5e-3 0 2e-2];
writeFirstDataFile(filename,N_goal,atom_types,simBox,atoms)

% define properties needed for writing HEAVY spheres data file 
filename = '../spheres_heavy.data';
rho = 100*1560; % particle density
sphereDensities = rho*ones(N_goal,1);
atoms = horzcat(id,atom_type,2*sphereRadii,sphereDensities,sphereCenters);
writeFirstDataFile(filename,N_goal,atom_types,simBox,atoms)

%% generate square punch
R_punch = box_x/20;
N_x = 24;
N_y = 24; 

box_x = 8e-3;
box_y = 8e-3;
box_z = 1e-2;
boxSize = [box_x, box_y, box_z]; 

% center point location of box
boxCenter = [0,0,box_z/2];

% gap width
w = 1.739738e-4;

[sphereCenters, sphereRadii] = generateSquarePunch(R_punch,boxSize,boxCenter,N_x,N_y,w);

% define properties needed for writing data file 
filename = '../sphere_punch.data';
N_punchMax = length(sphereRadii);
id = ((N_goal+1):1:(N_goal+N_punchMax))';
N_punch = length(id);
atom_type = 2*ones(N_punch,1);
rho = 1560; % particle density
sphereDensities = rho*ones(N_punch,1);
atoms = horzcat(id,atom_type,2*sphereRadii,sphereDensities,sphereCenters);
writeFirstDataFile(filename,N_punch,atom_types,simBox,atoms)

% tight punch
w = 1.739738e-5;

[sphereCenters, sphereRadii] = generateSquarePunch(R_punch,boxSize,boxCenter,N_x,N_y,w);

% define properties needed for writing data file 
filename = '../sphere_punch_tight.data';
N_punchMax = length(sphereRadii);
id = ((N_goal+1):1:(N_goal+N_punchMax))';
N_punch = length(id);
atom_type = 2*ones(N_punch,1);
rho = 1560; % particle density
sphereDensities = rho*ones(N_punch,1);
atoms = horzcat(id,atom_type,2*sphereRadii,sphereDensities,sphereCenters);
writeFirstDataFile(filename,N_punch,atom_types,simBox,atoms)


