clear all
close all
clc



force_disp = readmatrix('upperPunchDispForce.csv');


output_rate = round(1e-3/0.75e-7);
disp = force_disp(:,1);

n = length(disp);
new_length = floor(n / output_rate) * output_rate;

disp = disp(1:new_length);
disp = mean(reshape(disp,output_rate,[]),1);

create_moving_disk_stl_files(disp);