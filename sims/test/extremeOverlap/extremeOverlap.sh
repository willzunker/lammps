#!/bin/bash
set -euo pipefail

# path to lmp executable
LMP=../../../build/lmp

"${LMP}" < in.extremeOverlap_0
"${LMP}" < in.extremeOverlap_1
"${LMP}" < in.extremeOverlap_2
"${LMP}" < in.extremeOverlap_3
"${LMP}" < in.extremeOverlap_4
"${LMP}" < in.extremeOverlap_5

"${LMP}" < in.extremeOverlapFlipped_0
"${LMP}" < in.extremeOverlapFlipped_1
"${LMP}" < in.extremeOverlapFlipped_2
"${LMP}" < in.extremeOverlapFlipped_3
"${LMP}" < in.extremeOverlapFlipped_4
"${LMP}" < in.extremeOverlapFlipped_5
