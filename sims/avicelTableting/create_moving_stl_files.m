function create_moving_stl_files(input_stl, disp, stl_name)
    % Read the input STL file
    TR = stlread(input_stl);

    % Get vertices and faces from the triangulation object
    vertices = TR.Points;
    faces = TR.ConnectivityList;

    % Create STL files for each displacement
    for i = 1:length(disp)
        % Calculate the z-shift
        z_shift = disp(i);

        % Create a copy of the vertices and apply the z-shift
        shifted_vertices = vertices;
        shifted_vertices(:, 3) = shifted_vertices(:, 3) + z_shift;

        % Create a new triangulation object with shifted vertices
        shifted_TR = triangulation(faces, shifted_vertices);

        % Create the STL file
        stl_filename = sprintf('post/%s_%03d.stl', stl_name, i);
        stlwrite(shifted_TR, stl_filename);
    end
end