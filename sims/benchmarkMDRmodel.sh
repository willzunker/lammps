#!/bin/bash
set -euo pipefail

# path to lmp executable
LMP=../../../build/lmp

# non-adhesive simulations
cd uniaxialCompression
${LMP} < in.uniaxialCompression

cd ../dieCompaction
${LMP} < in.dieCompaction

cd ../triaxialCompaction
${LMP} < in.triaxialCompaction

# adhesive simulations
cd ../stickyParticleSandwich
${LMP} < in.stickyParticleSandwich

cd ../stickyUniaxialCompression
${LMP} < in.stickyUniaxialCompression
${LMP} < in.plasticStickyUniaxialCompression
