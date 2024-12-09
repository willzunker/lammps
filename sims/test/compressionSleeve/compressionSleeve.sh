#!/bin/bash
set -euo pipefail

# path to lmp executable
LMP=../../../build/lmp

# non-adhesive simulations
rm -f pairContactsTopCen.csv
rm -f pairContactsBotCen.csv
if [[ -d "post" ]]; then
    cd post
    rm -f *.vtk
    cd ..
fi
${LMP} < in.compressionSleeve
