% INPUT:
%   nWant:  Number of spheres
%   Width:  Dimension of 3d box as [1 x 3] double vector
%   Radius: Radius of spheres
% OUTPUT:
%   P:      [nWant x 3] matrix, centers

function [sphereCenters, sphereRadii] = placeCylSpheresWeibull(cylSize, cylCenter, N_goal, R_factor, R_min, R_max)
sphereCenters = zeros(N_goal, 3);
sphereRadii = zeros(N_goal,1);
N_cur = 0;
while  N_cur < N_goal
    radius = R_factor*(wblrnd(6, 0.69) + 4.5)*1e-6;
    if (radius >= R_min) && (radius <= R_max)
        sphereRadii(N_cur+1) = radius;
        N_cur = N_cur + 1
    end
end
place_r = cylSize(1) - max(sphereRadii);
place_z = cylSize(2) - 2*max(sphereRadii);
place_t = 2*pi;
i = 1;                  % Security break to avoid infinite loop
N_cur = 0;
while N_cur < N_goal && i < 3e6
    newSphereCenter_r = rand(1)*place_r;
    newSphereCenter_z = rand(1)*place_z + sphereRadii(N_cur+1);
    newSphereCenter_t = rand(1)*place_t;

    newSphereCenter_x = newSphereCenter_r*cos(newSphereCenter_t);
    newSphereCenter_y = newSphereCenter_r*sin(newSphereCenter_t);
    
    newSphereCenter = [newSphereCenter_x, newSphereCenter_y, newSphereCenter_z];
    if N_cur == 0
        N_cur = N_cur + 1;  % Append this point
        sphereCenters(N_cur, :) = newSphereCenter;
        continue
    end
    distance = sqrt(sum((sphereCenters(1:N_cur, :) - newSphereCenter) .^ 2, 2));
    touch = zeros(length(distance),1);
    for j = 1 : length(distance)
        if distance(j) <= (sphereRadii(N_cur+1) + sphereRadii(j))
            touch(j) = 1;
            break
        end
    end
    if sum(touch) == 0
        N_cur = N_cur + 1;  % Append this point
        sphereCenters(N_cur, :) = newSphereCenter;
    end

    i = i + 1;
    N_cur
end

packingFraction = sum(4/3*pi*sphereRadii.^3)./(pi*cylSize(1)^2*cylSize(2));
fprintf("The final packing fraction is: %f \n",packingFraction)

if N_cur < N_goal
    error('Cannot find wanted number of points in %d iterations.', i)
end

%plotSpheresCylinder(sphereCenters, sphereRadii)
end