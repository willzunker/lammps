#!/bin/bash

/Users/willzunker/lammps_mdr_develop/build/lmp < in.particleSandwich
mpirun --np 2 /Users/willzunker/lammps_mdr_develop/build/lmp -in in.particleSandwichParallel

/Users/willzunker/lammps_mdr_develop/build/lmp < in.particleSandwichCorner
mpirun --np 2 /Users/willzunker/lammps_mdr_develop/build/lmp -in in.particleSandwichCornerParallel
