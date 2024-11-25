% Plot spheres with specified centers and radii

function plotSpheresCylinder(sphereCenters, sphereRadii)
    figure()   
    hold on;
    [x, y, z] = sphere(); % Generate unit sphere coordinates
    
    for i = 1:size(sphereCenters, 1)
        % Scale the unit sphere to the desired radius and center it at the specified position
        xs = x * sphereRadii(i) + sphereCenters(i, 1);
        ys = y * sphereRadii(i) + sphereCenters(i, 2);
        zs = z * sphereRadii(i) + sphereCenters(i, 3);
        
        % Plot the scaled sphere
        surf(xs, ys, zs);
    end
    
    axis equal;
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    %xlim([(boxCenter(1)-boxSize(1)/2) (boxCenter(1)+boxSize(1)/2)])
    %ylim([(boxCenter(2)-boxSize(2)/2) (boxCenter(2)+boxSize(2)/2)])
    %zlim([(boxCenter(3)-boxSize(3)/2) (boxCenter(3)+boxSize(3)/2)])
    title('Spheres with Specific Radius');
    view(3);
    hold off;
end