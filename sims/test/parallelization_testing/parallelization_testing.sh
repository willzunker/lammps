#!/bin/bash

/Users/willzunker/lammps/build/lmp < in.E_Y_1
mpirun --np 1 /Users/willzunker/lammps/build/lmp -in in.E_Y_1_para1
mpirun --np 2 /Users/willzunker/lammps/build/lmp -in in.E_Y_1_para2
mpirun --np 4 /Users/willzunker/lammps/build/lmp -in in.E_Y_1_para4
mpirun --np 8 /Users/willzunker/lammps/build/lmp -in in.E_Y_1_para8