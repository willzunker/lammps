#!/bin/bash

# non-adhesive simulations
cd uniaxialCompression
/Users/willzunker/lammps/build/lmp < in.uniaxialCompression

cd ../dieCompaction
/Users/willzunker/lammps/build/lmp < in.dieCompaction

cd ../triaxialCompaction
/Users/willzunker/lammps/build/lmp < in.triaxialCompaction

# adhesive simulations
cd ../stickyParticleSandwich
 /Users/willzunker/lammps/build/lmp < in.stickyParticleSandwich

 cd ../stickyUniaxialCompression
 /Users/willzunker/lammps/build/lmp < in.stickyUniaxialCompression
 /Users/willzunker/lammps/build/lmp < in.plasticStickyUniaxialCompression