#!/bin/bash

# non-adhesive simulations
rm avgStresses.csv
cd post
rm *.vtk
cd ..
/Users/willzunker/lammps/build/lmp < in.avicelTableting