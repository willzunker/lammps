clear all
close all
clc

tic

%% initial parameters 
parameters.die_height = 0.01;
parameters.E = 1e6;        % young's modulus
parameters.nu = 0.3;       % poissons ratio
parameters.Y = 1e4;        % yield stress
parameters.gamma_p = 0.033;  % effective surface energy, particle-particle
parameters.gamma_w = 0.0;  % effective surface energy, particle-wall
parameters.psi_b = 0.5;    % bulk regime trigger
parameters.mu_s_p = 1.0;   % coefficent of friction, particle-particle
parameters.mu_s_w = 0.1;   % coefficent of friction, particle-wall
parameters.mu_r = 0.7;     % coefficent of rolling resistance
parameters.rho = 1560;     % density
parameters.Rmean = 0.33e-3; % mean radius
parameters.Rstd = parameters.Rmean/20;
parameters.dt_crit = sqrt(parameters.rho*parameters.Rmean^2/parameters.E);
parameters.dt = parameters.dt_crit*0.8;
parameters.skin_dis = 2.0*parameters.Rmean;
parameters.strain_target = 0.35;
parameters.powder_height = 0.6;
parameters.punch_stroke_prefactor = 0.6;

save('parameters.mat','parameters')

%% E/Y, gamma_w/(YRmean) variation for different Rstd/Rmean

gammapYRmean_o = parameters.gamma_p/(parameters.Rmean*parameters.Y);

 RstdRmean = [0,0.02,0.04,0.06,0.08,0.1];
% Rstd = parameters.Rmean*RstdRmean;
% 
% for k = 1 : length(Rstd)
% 
%     generateCylinderPacking(parameters.Rmean,Rstd(k))
% 
     EY = [10 20 40 80 160 320 640];
     gammapYRmean = [gammapYRmean_o/8, gammapYRmean_o/4, gammapYRmean_o/2, gammapYRmean_o, gammapYRmean_o*2, gammapYRmean_o*4, gammapYRmean_o*8];
% 
%     E = EY.*parameters.Y;
%     gammap = gammapYRmean.*(parameters.Y*parameters.Rmean);
     results_file = 'mech_prop';
%     for i = 1 : length(E)
%         parameters.dt_crit = sqrt(parameters.rho*parameters.Rmean^2/E(i));
%         parameters.dt = parameters.dt_crit*0.8;
%         for j = 1 : length(gammap)
%             for set_strain = 1 : 2
%                 if set_strain == 1
%                     % create script
%                     mech_prop_script(parameters,k,i,j,results_file,E(i),gammap(j));
% 
%                     % run script
%                     command = sprintf('cd %s && /Users/willzunker/lammps/build/lmp < in.%s_%d_%d%d', results_file, results_file,k,i,j);
%                     status_sim = system(command);
%                 else
%                     % correct settings to produce target axial strain
%                     force_disp_file = strcat('/Users/willzunker/lammps/sims/residual_stress_study/',results_file,'/',results_file,'_',num2str(k),'_',num2str(i),num2str(j),'_force_disp.csv');
%                     force_disp = abs(readmatrix(force_disp_file));
%                     max_force = max(force_disp(:,2));
%                     threshold = 0.03*max_force;
%                     indices = find(force_disp(:,2) > threshold);
%                     first_index = indices(1);
%                     disp_touch = force_disp(first_index,1);
%                     parameters.powder_height = parameters.die_height-disp_touch;
%                     parameters.punch_stroke_prefactor = (parameters.die_height-parameters.powder_height*(1-parameters.strain_target))/parameters.die_height;
% 
%                     % create script with correct axial strain
%                     mech_prop_script(parameters,k,i,j,results_file,E(i),gammap(j));
% 
%                     % run script
%                     command = sprintf('cd %s && /Users/willzunker/lammps/build/lmp < in.%s_%d_%d%d', results_file, results_file,k,i,j);
%                     status_sim = system(command);
% 
%                     % plot stress-strain
%                     %plot_stress_strain(parameters,k,i,j,results_file,parameters.Y);
%                 end
%             end
%         end
%     end
% end
% 
 plot_residual_stresses_mech_prop(results_file,parameters,EY,gammapYRmean,RstdRmean)


%% Wall particle friction variation

generateCylinderPacking(parameters.Rmean,parameters.Rstd)

mu_s_w = linspace(0,1,50);
%mu_s_w = 0.1;

results_file = 'mu_s_w';
for i = 1 : length(mu_s_w)
    for set_strain = 1 : 2
        if set_strain == 1
            % create script
            mu_s_w_script(parameters,i,results_file,mu_s_w(i));

            % run script
            command = sprintf('cd %s && /Users/willzunker/lammps/build/lmp < in.%s_%d', results_file, results_file, i);
            status_sim = system(command);
        else
            % correct settings to produce target axial strain
            force_disp_file = strcat('/Users/willzunker/lammps/sims/residual_stress_study/',results_file,'/',results_file,'_',num2str(i),'_force_disp.csv'); 
            force_disp = abs(readmatrix(force_disp_file));
            max_force = max(force_disp(:,2));
            threshold = 0.03*max_force;
            indices = find(force_disp(:,2) > threshold);
            first_index = indices(1);
            disp_touch = force_disp(first_index,1);
            parameters.powder_height = parameters.die_height-disp_touch;
            parameters.punch_stroke_prefactor = (parameters.die_height-parameters.powder_height*(1-parameters.strain_target))/parameters.die_height; 

            % create script with correct axial strain
            mu_s_w_script(parameters,i,results_file,mu_s_w(i));

            % run script
            command = sprintf('cd %s && /Users/willzunker/lammps/build/lmp < in.%s_%d', results_file, results_file, i);
            status_sim = system(command);

            % plot stress-strain
            %plot_stress_strain(parameters,0,i,0,results_file,parameters.Y);
        end
    end
end

plot_residual_stresses_mu_s_w(results_file,parameters.Y,mu_s_w)

%% Particle-particle friction variation 

mu_s_p = linspace(0,1,50);

results_file = 'mu_s_p';
for i = 1 : length(mu_s_p)
    for set_strain = 1 : 2
        if set_strain == 1
            % create script
            mu_s_p_script(parameters,i,results_file,mu_s_p(i));

            % run script
            command = sprintf('cd %s && /Users/willzunker/lammps/build/lmp < in.%s_%d', results_file, results_file, i);
            status_sim = system(command);
        else
            % correct settings to produce target axial strain
            force_disp_file = strcat('/Users/willzunker/lammps/sims/residual_stress_study/',results_file,'/',results_file,'_',num2str(i),'_force_disp.csv'); 
            force_disp = abs(readmatrix(force_disp_file));
            max_force = max(force_disp(:,2));
            threshold = 0.03*max_force;
            indices = find(force_disp(:,2) > threshold);
            first_index = indices(1);
            disp_touch = force_disp(first_index,1);
            parameters.powder_height = parameters.die_height-disp_touch;
            parameters.punch_stroke_prefactor = (parameters.die_height-parameters.powder_height*(1-parameters.strain_target))/parameters.die_height; 

            % create script with correct axial strain
            mu_s_p_script(parameters,i,results_file,mu_s_p(i));

            % run script
            command = sprintf('cd %s && /Users/willzunker/lammps/build/lmp < in.%s_%d', results_file, results_file, i);
            status_sim = system(command);

            % plot stress-strain
            %plot_stress_strain(parameters,0,i,0,results_file,parameters.Y);
        end
    end
end

plot_residual_stresses_mu_s_p(results_file,parameters.Y,mu_s_p)

%% Rolling resistance coefficent variation 

mu_r = linspace(0,1,50);

results_file = 'mu_r';
for i = 1 : length(mu_r)
    for set_strain = 1 : 2
        if set_strain == 1
            % create script
            mu_r_script(parameters,i,results_file,mu_r(i));

            % run script
            command = sprintf('cd %s && /Users/willzunker/lammps/build/lmp < in.%s_%d', results_file, results_file, i);
            status_sim = system(command);
        else
            % correct settings to produce target axial strain
            force_disp_file = strcat('/Users/willzunker/lammps/sims/residual_stress_study/',results_file,'/',results_file,'_',num2str(i),'_force_disp.csv'); 
            force_disp = abs(readmatrix(force_disp_file));
            max_force = max(force_disp(:,2));
            threshold = 0.03*max_force;
            indices = find(force_disp(:,2) > threshold);
            first_index = indices(1);
            disp_touch = force_disp(first_index,1);
            parameters.powder_height = parameters.die_height-disp_touch;
            parameters.punch_stroke_prefactor = (parameters.die_height-parameters.powder_height*(1-parameters.strain_target))/parameters.die_height; 

            % create script with correct axial strain
            mu_r_script(parameters,i,results_file,mu_r(i));

            % run script
            command = sprintf('cd %s && /Users/willzunker/lammps/build/lmp < in.%s_%d', results_file, results_file, i);
            status_sim = system(command);

            % plot stress-strain
            %plot_stress_strain(parameters,0,i,0,results_file,parameters.Y);
        end
    end
end

plot_residual_stresses_mu_r(results_file,parameters.Y,mu_r)

%% change final strain level

strain = linspace(0.05,0.4,50);

results_file = 'strain';
for i = 1 : length(strain)
    for set_strain = 1 : 2
        if set_strain == 1
            % create script
            strain_script(parameters,i,results_file);

            % run script
            command = sprintf('cd %s && /Users/willzunker/lammps/build/lmp < in.%s_%d', results_file, results_file, i);
            status_sim = system(command);
        else
            % correct settings to produce target axial strain
            force_disp_file = strcat('/Users/willzunker/lammps/sims/residual_stress_study/',results_file,'/',results_file,'_',num2str(i),'_force_disp.csv'); 
            force_disp = abs(readmatrix(force_disp_file));
            max_force = max(force_disp(:,2));
            threshold = 0.03*max_force;
            indices = find(force_disp(:,2) > threshold);
            first_index = indices(1);
            disp_touch = force_disp(first_index,1);
            parameters.powder_height = parameters.die_height-disp_touch;
            parameters.punch_stroke_prefactor = (parameters.die_height-parameters.powder_height*(1-strain(i)))/parameters.die_height; 

            % create script with correct axial strain
            strain_script(parameters,i,results_file);

            % run script
            command = sprintf('cd %s && /Users/willzunker/lammps/build/lmp < in.%s_%d', results_file, results_file, i);
            status_sim = system(command);

            % plot stress-strain
            %plot_stress_strain(parameters,0,i,0,results_file,parameters.Y);
        end
    end
end

plot_residual_stresses_strain(results_file,parameters.Y,strain)

%% change psi_b

psi_b = linspace(0.4,0.8,50);

results_file = 'psi_b';
for i = 1 : length(psi_b)
    for set_strain = 1 : 2
        if set_strain == 1
            % create script
            psi_b_script(parameters,i,results_file,psi_b(i));

            % run script
            command = sprintf('cd %s && /Users/willzunker/lammps/build/lmp < in.%s_%d', results_file, results_file, i);
            status_sim = system(command);
        else
            % correct settings to produce target axial strain
            force_disp_file = strcat('/Users/willzunker/lammps/sims/residual_stress_study/',results_file,'/',results_file,'_',num2str(i),'_force_disp.csv'); 
            force_disp = abs(readmatrix(force_disp_file));
            max_force = max(force_disp(:,2));
            threshold = 0.03*max_force;
            indices = find(force_disp(:,2) > threshold);
            first_index = indices(1);
            disp_touch = force_disp(first_index,1);
            parameters.powder_height = parameters.die_height-disp_touch;
            parameters.punch_stroke_prefactor = (parameters.die_height-parameters.powder_height*(1-parameters.strain_target))/parameters.die_height; 

            % create script with correct axial strain
            psi_b_script(parameters,i,results_file,psi_b(i));

            % run script
            command = sprintf('cd %s && /Users/willzunker/lammps/build/lmp < in.%s_%d', results_file, results_file, i);
            status_sim = system(command);

            % plot stress-strain
            %plot_stress_strain(parameters,0,i,0,results_file,parameters.Y);
        end
    end
end

plot_residual_stresses_psi_b(results_file,parameters.Y,psi_b)

toc