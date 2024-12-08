#!/bin/bash

# non-adhesive simulations
rm pairContactsTopCen.csv
rm pairContactsBotCen.csv
cd post
rm *.vtk
cd ..
/Users/willzunker/lammps/build/lmp < in.compressionSleeveFree