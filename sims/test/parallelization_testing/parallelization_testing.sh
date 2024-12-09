#!/bin/bash
set -euo pipefail

# path to lmp executable
LMP=../../../build/lmp

"${LMP}" < in.E_Y_1
mpirun --np 1 "${LMP}" -in in.E_Y_1_para1
mpirun --np 2 "${LMP}" -in in.E_Y_1_para2
mpirun --np 4 "${LMP}" -in in.E_Y_1_para4
