clear all
close all
clc

    % Parameters of the box
    center = [0, 0, 5e-3];
    width = 8e-3;
    height = 1e-2;
    
    % Calculate box dimensions
    half_width = width / 2;
    half_height = height / 2;

    % Define the 8 vertices of the box
    vertices = [
        center(1) - half_width, center(2) - half_width, center(3) - half_height;
        center(1) + half_width, center(2) - half_width, center(3) - half_height;
        center(1) + half_width, center(2) + half_width, center(3) - half_height;
        center(1) - half_width, center(2) + half_width, center(3) - half_height;
        center(1) - half_width, center(2) - half_width, center(3) + half_height;
        center(1) + half_width, center(2) - half_width, center(3) + half_height;
        center(1) + half_width, center(2) + half_width, center(3) + half_height;
        center(1) - half_width, center(2) + half_width, center(3) + half_height;
    ];

    % Define the 10 faces of the box (without the top face)
    faces = [
        1 2 3; 1 3 4; % Bottom face
        1 2 6; 1 6 5; % Side face 1
        2 3 7; 2 7 6; % Side face 2
        3 4 8; 3 8 7; % Side face 3
        4 1 5; 4 5 8; % Side face 4
        5 6 7; 5 7 8; % Top face (inside, optional)
    ];

    % Open a new STL file
    stl_file = fopen('confining_box.stl', 'w');
    fprintf(stl_file, 'solid box\n');
    
    % Write each face to the STL file
    for i = 1:size(faces, 1)
        v1 = vertices(faces(i, 1), :);
        v2 = vertices(faces(i, 2), :);
        v3 = vertices(faces(i, 3), :);
        fprintf(stl_file, '  facet normal 0 0 0\n');
        fprintf(stl_file, '    outer loop\n');
        fprintf(stl_file, '      vertex %.6f %.6f %.6f\n', v1);
        fprintf(stl_file, '      vertex %.6f %.6f %.6f\n', v2);
        fprintf(stl_file, '      vertex %.6f %.6f %.6f\n', v3);
        fprintf(stl_file, '    endloop\n');
        fprintf(stl_file, '  endfacet\n');
    end
    
    fprintf(stl_file, 'endsolid box\n');
    
    % Close the STL file
    fclose(stl_file);
