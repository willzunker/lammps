#!/bin/bash

/Users/willzunker/lammps/build/lmp < in.twoParticle_1
/Users/willzunker/lammps/build/lmp < in.twoParticle_0_75
/Users/willzunker/lammps/build/lmp < in.twoParticle_0_5
/Users/willzunker/lammps/build/lmp < in.twoParticle_0_25
/Users/willzunker/lammps/build/lmp < in.twoParticle_0

/Users/willzunker/lammps/build/lmp < in.twoParticleFlipped_1
/Users/willzunker/lammps/build/lmp < in.twoParticleFlipped_0_75
/Users/willzunker/lammps/build/lmp < in.twoParticleFlipped_0_5
/Users/willzunker/lammps/build/lmp < in.twoParticleFlipped_0_25
/Users/willzunker/lammps/build/lmp < in.twoParticleFlipped_0