/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef ANGLE_CLASS
// clang-format off
AngleStyle(amoeba,AngleAmoeba);
// clang-format on
#else

#ifndef LMP_ANGLE_AMOEBA_H
#define LMP_ANGLE_AMOEBA_H

#include "angle.h"

namespace LAMMPS_NS {

class AngleAmoeba : public Angle {
 public:
  AngleAmoeba(class LAMMPS *);
  virtual ~AngleAmoeba();
  void compute(int, int);
  void coeff(int, char **);
  double equilibrium_angle(int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_data(FILE *);
  double single(int, int, int, int);

 protected:
  int *pflag;
  double *theta0, *k2, *k3, *k4, *k5, *k6;

  void anglep(int, int, int, int, int);
  void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Incorrect args for angle coefficients

Self-explanatory.  Check the input script or data file.

*/
