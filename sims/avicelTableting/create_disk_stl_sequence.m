clear all
close all
clc



disp = readmatrix('punch_disp_STL_generation.csv');

dt = 0.3e-7;
output_rate = round(1e-4/dt);
disp_upper = disp(:,1);
disp_lower = disp(:,2);

n = length(disp);
new_length = floor(n / output_rate) * output_rate;

disp_upper = disp_upper(1:new_length);
disp_lower = disp_lower(1:new_length);

disp_upper = mean(reshape(disp_upper,output_rate,[]),1);
disp_lower = mean(reshape(disp_lower,output_rate,[]),1);

zo_upper = 0.01;
zo_lower = 0.0;

%create_moving_disk_stl_files(zo_upper,disp_upper,'disk_upper');
%create_moving_disk_stl_files(zo_lower,disp_lower,'disk_lower');

create_moving_stl_files('upper_punch.stl', disp_upper, 'upper_punch');
create_moving_stl_files('lower_punch.stl', disp_lower, 'lower_punch');