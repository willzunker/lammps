function writeDataFile(filename, N, atoms)
    % Open file for writing
    fid = fopen(filename, 'w');
    
    % Write header
    fprintf(fid, '#LAMMPS data file created by matlab.\n');
    fprintf(fid, '%d atoms\n\n', N);
    
    % Write atoms
    fprintf(fid, 'Atoms\n\n');
    for i = 1:N
        fprintf(fid, '%d %d %f %f %f %f %f\n', atoms(i,:));
    end
    
    % Close file
    fclose(fid);
end
