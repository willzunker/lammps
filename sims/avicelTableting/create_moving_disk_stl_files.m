function create_moving_disk_stl_files(disp)
    % Parameters
    radius = 4e-3;   % Radius of the disk
    initial_z = 0.01; % Initial z-coordinate of the disk center
    num_points = 100; % Number of points to define the disk's circumference

    % Generate the points for the disk
    theta = linspace(0, 2*pi, num_points);
    x = radius * cos(theta);
    y = radius * sin(theta);
    z = zeros(1, num_points);
    
    % Include the center of the disk as the first vertex
    vertices = [0, 0, 0; x', y', z'];
    
    % Define the faces (triangles) of the disk
    faces = [(2:num_points)', (3:num_points+1)', ones(num_points-1, 1)];
    faces = [faces; num_points+1, 2, 1]; % Close the disk

    % Create STL files for each displacement
    for i = 1:length(disp)
        % Current z-coordinate of the disk center
        current_z = initial_z + disp(i);

        % Update z-coordinates of the vertices
        vertices(:, 3) = [0; z'] + current_z;

        % Create the STL file
        stl_filename = sprintf('post/disk_%03d.stl', i);
        create_stl_file(stl_filename, vertices, faces);
    end
end

function create_stl_file(filename, vertices, faces)
    % Create an STL file with the given vertices and faces
    num_faces = size(faces, 1);
    fid = fopen(filename, 'w');
    fprintf(fid, 'solid disk\n');
    
    for i = 1:num_faces
        v1 = vertices(faces(i, 1), :);
        v2 = vertices(faces(i, 2), :);
        v3 = vertices(faces(i, 3), :);
        normal = cross(v2 - v1, v3 - v1);
        normal = normal / norm(normal);
        fprintf(fid, 'facet normal %e %e %e\n', normal);
        fprintf(fid, 'outer loop\n');
        fprintf(fid, 'vertex %e %e %e\n', v1);
        fprintf(fid, 'vertex %e %e %e\n', v2);
        fprintf(fid, 'vertex %e %e %e\n', v3);
        fprintf(fid, 'endloop\n');
        fprintf(fid, 'endfacet\n');
    end
    
    fprintf(fid, 'endsolid disk\n');
    fclose(fid);
end

