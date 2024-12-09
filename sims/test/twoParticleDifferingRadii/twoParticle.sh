#!/bin/bash
set -euo pipefail

# path to lmp executable
LMP=../../../build/lmp

"${LMP}" < in.twoParticle_1
"${LMP}" < in.twoParticle_0_75
"${LMP}" < in.twoParticle_0_5
"${LMP}" < in.twoParticle_0_25
"${LMP}" < in.twoParticle_0

"${LMP}" < in.twoParticleFlipped_1
"${LMP}" < in.twoParticleFlipped_0_75
"${LMP}" < in.twoParticleFlipped_0_5
"${LMP}" < in.twoParticleFlipped_0_25
"${LMP}" < in.twoParticleFlipped_0
