#!/bin/bash
set -euo pipefail

# path to lmp executable
LMP=../../build/lmp

# serial vs parallel
${LMP} < in.avicelTableting
#mpirun --np 4 ${LMP} -in in.avicelTableting
