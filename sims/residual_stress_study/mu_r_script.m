function mu_r_script(parameters,i,results_file,mu_r)

ID = num2str(i);
filename = strcat(results_file,'/in.mu_r_',ID);
dump_name = strcat('mu_r_',ID,'.dump');
force_disp_name = strcat('mu_r_',ID,'_force_disp.csv');
avg_stresses_name = strcat('mu_r_',ID,'_avg_stresses.csv');

% Open the file for writing
fileID = fopen(filename, 'w');

% Write the LAMMPS input script to the file
fprintf(fileID, '############################### SIMULATION SETTINGS ###################################################\n\n');
fprintf(fileID, 'atom_style sphere \n');
fprintf(fileID, 'atom_modify map array \n');
fprintf(fileID, 'comm_modify vel yes\n');
fprintf(fileID, 'units si\n');
fprintf(fileID, 'newton off\n');
fprintf(fileID, 'neighbor %e bin\n',parameters.skin_dis);
fprintf(fileID, 'neigh_modify every 1 delay 0 check yes \n');
fprintf(fileID, 'timestep %e\n\n',parameters.dt);

fprintf(fileID, '######################### SIMULATION BOUNDING BOX, INTEGRATION, AND, GRAVITY ###########################\n\n');
fprintf(fileID, 'boundary f f f\n');
fprintf(fileID, 'read_data ../spheres.data\n\n');

fprintf(fileID, '######################################### ADD DIE AND ATOM PARAMETERIZATION ##############################################\n\n');
fprintf(fileID, 'variable atomRadius equal %e\n',parameters.Rmean);
fprintf(fileID, 'variable atomDiameter equal 2*${atomRadius}\n');
fprintf(fileID, 'variable atomDensity equal 1560 \n');
fprintf(fileID, 'variable dieRadius equal 4e-3\n');
fprintf(fileID, 'variable dieHeight equal 1e-2\n\n');

fprintf(fileID, '########################### PARTICLE MATERIAL PROPERTIES AND FORCE MODEL ###############################\n\n');
fprintf(fileID, 'pair_style granular\n');
fprintf(fileID, '# mdr = E, nu, Y, gamma, psi_b, CoR\n');
fprintf(fileID, '# linear_history = k_t, x_gamma,t, mu_s \n');
fprintf(fileID, 'variable YoungsModulus equal %e\n', parameters.E);
fprintf(fileID, 'variable PoissonsRatio equal %f\n', parameters.nu);
fprintf(fileID, 'variable YieldStress equal %e\n',parameters.Y);
fprintf(fileID, 'variable SurfaceEnergy equal %e\n',parameters.gamma_p);
fprintf(fileID, 'variable SurfaceEnergyWall equal %e\n',parameters.gamma_w);
fprintf(fileID, 'variable psi_b equal %f\n',parameters.psi_b);
fprintf(fileID, 'variable CoR equal 0.5\n');
fprintf(fileID, 'variable kt equal 2/7*${YoungsModulus}*${atomRadius}\n');
fprintf(fileID, 'variable kt_wall equal 2/7*${YoungsModulus}*${atomRadius}\n');
fprintf(fileID, 'variable xgammat equal 0.0\n');
fprintf(fileID, 'variable mu_s equal %f\n',parameters.mu_s_p);
fprintf(fileID, 'variable mu_s_wall equal %f\n',parameters.mu_s_w);
fprintf(fileID, 'variable mu_roll equal %f\n',mu_r);
fprintf(fileID, 'variable k_roll equal 2.25*${mu_roll}*${mu_roll}*${YoungsModulus}*${atomRadius}\n');
fprintf(fileID, 'variable gamma_roll equal 0.0\n\n');

fprintf(fileID, 'pair_coeff * * mdr ${YoungsModulus} ${PoissonsRatio} ${YieldStress} ${SurfaceEnergy} ${psi_b} ${CoR} damping none tangential linear_history ${kt} ${xgammat} ${mu_s} rolling sds ${k_roll} ${gamma_roll} ${mu_roll}\n');
fprintf(fileID, '#pair_coeff * * mdr ${YoungsModulus} ${PoissonsRatio} ${YieldStress} ${SurfaceEnergy} ${psi_b} ${CoR} tangential linear_history ${kt} ${xgammat} ${mu_s} damping none\n');
fprintf(fileID, '#pair_coeff * * mdr ${YoungsModulus} ${PoissonsRatio} ${YieldStress} ${SurfaceEnergy} ${psi_b} ${CoR} tangential linear_nohistory 0.0 0.0 damping none\n\n');

fprintf(fileID, '######################################### ADD DIE AND PUNCH WALLS ################################################\n\n');
fprintf(fileID, 'variable disp_upper equal 0.0\n');
fprintf(fileID, 'variable disp_lower equal 0.0\n\n');

fprintf(fileID, 'variable wall_contact_string string "granular mdr ${YoungsModulus} ${PoissonsRatio} ${YieldStress} ${SurfaceEnergyWall} ${psi_b} ${CoR} damping none tangential linear_history ${kt_wall} ${xgammat} ${mu_s_wall} rolling sds ${k_roll} ${gamma_roll} ${mu_roll}"\n');
fprintf(fileID, '#variable wall_contact_string string "granular mdr ${YoungsModulus} ${PoissonsRatio} ${YieldStress} ${SurfaceEnergyWall} ${psi_b} ${CoR} tangential linear_history ${kt} ${xgammat} ${mu_s} damping none"\n');
fprintf(fileID, '#variable wall_contact_string string "granular mdr ${YoungsModulus} ${PoissonsRatio} ${YieldStress} ${SurfaceEnergyWall} ${psi_b} ${CoR} tangential linear_nohistory 0.0 0.0 damping none"\n\n');

fprintf(fileID, 'variable dieHeight2 equal 2*${dieHeight}\n\n');

fprintf(fileID, 'region lowerPunch plane 0 0 0 0 0 1 side in units box move NULL NULL v_disp_lower units box\n');
fprintf(fileID, 'region upperPunch plane 0 0 ${dieHeight} 0 0 -1 side in move NULL NULL v_disp_upper units box\n');
fprintf(fileID, 'region die cylinder z 0 0 ${dieRadius} 0 ${dieHeight2} side in units box\n\n');

fprintf(fileID, 'fix lowerPunch all wall/gran/region ${wall_contact_string} region lowerPunch contacts\n');
fprintf(fileID, 'fix upperPunch all wall/gran/region ${wall_contact_string} region upperPunch contacts\n');
fprintf(fileID, 'fix die all wall/gran/region ${wall_contact_string} region die contacts \n\n');

fprintf(fileID, '##################################### INSERT PARTICLES ####################################################\n\n');
fprintf(fileID, 'fix 1 all nve/sphere\n');
fprintf(fileID, 'fix grav all gravity 9.81 vector 0 0 -1\n\n');

fprintf(fileID, '######################################## SCREEN OUTPUT  ####################################################\n\n');
fprintf(fileID, 'compute 1 all erotate/sphere\n');
fprintf(fileID, 'thermo_style custom dt step atoms ke vol\n');
fprintf(fileID, 'thermo 100\n');
fprintf(fileID, 'thermo_modify lost ignore norm no\n\n');

fprintf(fileID, '##################################### SET UP DUMP OUTPUTS  ####################################################\n\n');
fprintf(fileID, 'compute ke all ke/atom\n');
fprintf(fileID, 'variable output_rate equal round(1e-3/dt)\n');
fprintf(fileID, 'dump dumpParticles all custom ${output_rate} %s id type mass diameter x y z vx vy vz fx fy fz c_ke\n',dump_name);

fprintf(fileID, '#variable az_upperPunch atom f_upperPunch[7]\n');
fprintf(fileID, '#variable afz_upperPunch atom f_upperPunch[4]\n');

fprintf(fileID, 'run 0\n\n');

fprintf(fileID, 'compute sigmaxx all property/atom d_sigmaxx\n');
fprintf(fileID, 'compute sigmayy all property/atom d_sigmayy\n');
fprintf(fileID, 'compute sigmazz all property/atom d_sigmazz\n');
fprintf(fileID, 'compute Velas all property/atom d_Velas\n\n');

fprintf(fileID, 'compute sigmaxx_ave all reduce ave c_sigmaxx\n');
fprintf(fileID, 'compute sigmayy_ave all reduce ave c_sigmayy\n');
fprintf(fileID, 'compute sigmazz_ave all reduce ave c_sigmazz\n');
fprintf(fileID, 'compute Velas_sum all reduce sum c_Velas\n\n');

fprintf(fileID, 'variable sxx_ave equal c_sigmaxx_ave\n');
fprintf(fileID, 'variable syy_ave equal c_sigmayy_ave\n');
fprintf(fileID, 'variable szz_ave equal c_sigmazz_ave\n');
fprintf(fileID, 'variable Vparticles equal c_Velas_sum\n\n');

fprintf(fileID, 'fix log all print 1 "${sxx_ave} ${syy_ave} ${szz_ave} ${Vparticles}" file %s screen no\n\n',avg_stresses_name);

fprintf(fileID, 'compute avgPunchForce all reduce sum f_upperPunch[4]\n');
fprintf(fileID, 'variable avgPunchForce equal c_avgPunchForce\n\n');

fprintf(fileID, 'fix print1 all print 1 "${disp_upper} ${avgPunchForce}" file %s screen no \n\n',force_disp_name);

fprintf(fileID, '######################################### RUN SIMULATION ##########################################\n\n');

fprintf(fileID, 'variable upper_punch_stroke equal %f*${dieHeight}\n',parameters.punch_stroke_prefactor);
fprintf(fileID, 'variable vel_upper equal 0.25\n\n');

fprintf(fileID, 'variable settling_steps equal round(0.02/dt)\n');
fprintf(fileID, 'variable compression_steps equal 2*round(${upper_punch_stroke}/${vel_upper}/dt)\n');
fprintf(fileID, 'variable ejection_steps equal ${compression_steps}\n');
fprintf(fileID, 'variable free_float_steps equal round(0.025/dt)\n');
fprintf(fileID, 'variable total_steps equal ${settling_steps}+${compression_steps}+${ejection_steps}+${free_float_steps}\n\n');

fprintf(fileID, 'print "Total steps = ${total_steps}"\n\n');

fprintf(fileID, '##### SETTLING #####\n\n');
fprintf(fileID, 'run ${settling_steps}\n\n');

fprintf(fileID, '##### Compression & Release #####\n\n');
fprintf(fileID, 'variable punch_frequency equal PI/2/(dt*${compression_steps}/2)\n');
fprintf(fileID, 'variable disp_upper equal -${upper_punch_stroke}*sin(${punch_frequency}*elapsed*dt)\n');
fprintf(fileID, 'variable short_release equal round(${compression_steps}*1.0)\n');
fprintf(fileID, 'run ${short_release}\n\n');

%fprintf(fileID, '##### EJECTION #####\n\n');
%fprintf(fileID, 'variable punch_frequency equal PI/2/(dt*${ejection_steps})\n');
%fprintf(fileID, 'variable disp_lower equal ${dieHeight}*sin(${punch_frequency}*elapsed*dt)\n');
%fprintf(fileID, 'variable disp_upper equal 0.9*v_disp_lower\n');
%fprintf(fileID, 'run ${ejection_steps}\n\n');
%
%fprintf(fileID, '##### FREE FLOAT #####\n\n');
%fprintf(fileID, 'variable disp_lower equal ${dieHeight}\n');
%fprintf(fileID, 'variable disp_upper equal ${dieHeight}*0.9\n');
%fprintf(fileID, 'run ${free_float_steps}\n');

% Close the file
fclose(fileID);

end