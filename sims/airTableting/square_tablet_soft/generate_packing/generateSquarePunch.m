function [sphereCenters, sphereRadii] = generateSquarePunch(R_punch,boxSize,boxCenter,N_x,N_y,w)
box_x = boxSize(1);
box_y = boxSize(2);
box_z = boxSize(3);

z_punch = boxCenter(3)+box_z/2;

N = (N_x-1)*(N_y-1);

sphereCenters = zeros(N,3);
sphereRadii = R_punch*ones(N,1);

dx = box_x/N_x;
dy = box_y/N_y;

z_c = z_punch + R_punch;

k = 1;
for i = 1 : N_x-1
    for j = 1 : N_y-1
        if i == 1 && j == 1
            sphereCenters(k,1:3) = [i*dx+w,j*dy+w,z_c];
        elseif i == 1 && j == N_y-1
            sphereCenters(k,1:3) = [i*dx+w,j*dy-w,z_c];
        elseif i == N_x-1 && j == 1
            sphereCenters(k,1:3) = [i*dx-w,j*dy+w,z_c];
        elseif i == N_x-1 && j == N_y-1
            sphereCenters(k,1:3) = [i*dx-w,j*dy-w,z_c];
        elseif i == 1
            sphereCenters(k,1:3) = [i*dx+w,j*dy,z_c];
        elseif j == 1
            sphereCenters(k,1:3) = [i*dx,j*dy+w,z_c];
        elseif i == N_x-1
            sphereCenters(k,1:3) = [i*dx-w,j*dy,z_c];
        elseif j == N_y-1
            sphereCenters(k,1:3) = [i*dx,j*dy-w,z_c];
        else
            sphereCenters(k,1:3) = [i*dx,j*dy,z_c];
        end
        k = k+1;
    end
end

boxTranslate = boxCenter - boxSize/2;
sphereCenters = sphereCenters + boxTranslate;

plotSpheres(sphereCenters, sphereRadii, boxSize, boxCenter)

end

