#!/bin/bash
set -euo pipefail

# path to lmp executable
LMP=../../../build/lmp

cd test

# non-adhesive simulations
cd uniaxialCompression
mkdir -p post
${LMP} < in.uniaxialCompression

cd ../dieCompaction
mkdir -p post
${LMP} < in.dieCompaction

cd ../triaxialCompaction
mkdir -p post
${LMP} < in.triaxialCompaction

# adhesive simulations
cd ../stickyParticleSandwich
mkdir -p post
${LMP} < in.stickyParticleSandwich

cd ../stickyUniaxialCompression
mkdir -p post
${LMP} < in.stickyUniaxialCompression
${LMP} < in.plasticStickyUniaxialCompression
