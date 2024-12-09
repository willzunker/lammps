#!/bin/bash
set -euo pipefail

# path to lmp executable
LMP=../../../build/lmp

"${LMP}" < in.particleSandwich
mpirun --np 2 "${LMP}" -in in.particleSandwichParallel

"${LMP}" < in.particleSandwichCorner
mpirun --np 2 "${LMP}" -in in.particleSandwichCornerParallel
