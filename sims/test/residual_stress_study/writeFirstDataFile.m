function writeFirstDataFile(filename, N, atom_types, bounds, atoms)
    % Open file for writing
    fid = fopen(filename, 'w');
    
    % Write header
    fprintf(fid, '#LAMMPS data file created by matlab.\n');
    fprintf(fid, '%d atoms\n\n', N);
    fprintf(fid, '%d atom types\n\n', atom_types);
    
    % Write box bounds
    fprintf(fid, '%f %f xlo xhi\n', bounds(1), bounds(2));
    fprintf(fid, '%f %f ylo yhi\n', bounds(3), bounds(4));
    fprintf(fid, '%f %f zlo zhi\n\n', bounds(5), bounds(6));
    
    % Write atoms
    fprintf(fid, 'Atoms\n\n');
    for i = 1:N
        fprintf(fid, '%d %d %f %f %f %f %f\n', atoms(i,:));
    end
    
    % Close file
    fclose(fid);
end
