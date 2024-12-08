#!/bin/bash

# non-adhesive simulations
rm avgStresses.csv
/Users/willzunker/lammps/build/lmp < in.avicelTableting
#mpirun --np 9 /Users/willzunker/lammps/build/lmp -in in.avicelTableting