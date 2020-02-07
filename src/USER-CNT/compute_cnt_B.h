/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Contributing author: Maxim Shugaev (UVA), mvs9t@virginia.edu
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(cnt/B,ComputeCNT_B)

#else

#ifndef LMP_COMPUTE_CNT_B_ATOM_H
#define LMP_COMPUTE_CNT_B_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeCNT_B : public Compute {
 public:
  ComputeCNT_B(class LAMMPS *, int, char **);
  ~ComputeCNT_B();
  void init() {}
  void compute_peratom();
  double memory_usage();

 private:
  int nmax;
  double *buckling;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal compute cnt/B command

*/