# preset that turns on a wide range of packages, some of which require
# external libraries. Compared to all_on.cmake some more unusual packages
# are removed. The resulting binary should be able to run most inputs.

set(ALL_PACKAGES
  ASPHERE
  BOCS
  BODY
  BROWNIAN
  CG-DNA
  CG-SDK
  CLASS2
  COLLOID
  COLVARS
  COMPRESS
  CORESHELL
  DIELECTRIC
  DIFFRACTION
  DIPOLE
  DPD-BASIC
  DPD-MESO
  DPD-REACT
  DPD-SMOOTH
  DRUDE
  EFF
  FEP
  GRANULAR
  KSPACE
  MACHDYN
  MANYBODY
  MC
  MEAM
  MISC
  ML-IAP
  ML-SNAP
  MOFFF
  MOLECULE
  OPENMP
  OPT
  PERI
  PHONON
  PLUGIN
  POEMS
  PYTHON
  QEQ
  REACTION
  REAXFF
  REPLICA
  RIGID
  SHOCK
  SPH
  SPIN
  SRD
  UEF
  USER-MISC
  VORONOI
  YAFF)

foreach(PKG ${ALL_PACKAGES})
  set(PKG_${PKG} ON CACHE BOOL "" FORCE)
endforeach()

set(BUILD_TOOLS ON CACHE BOOL "" FORCE)
