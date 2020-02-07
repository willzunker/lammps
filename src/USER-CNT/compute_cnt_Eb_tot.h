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

ComputeStyle(cnt/Eb_tot,ComputeCNT_Eb_tot)

#else

#ifndef LMP_COMPUTE_CNT_EB_H
#define LMP_COMPUTE_CNT_EB_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeCNT_Eb_tot : public Compute {
 public:
  ComputeCNT_Eb_tot(class LAMMPS *, int, char **);
  ~ComputeCNT_Eb_tot() {}
  void init() {}
  double compute_scalar();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal compute cnt/Eb_tot command

Incorrect argument list in the compute init.

E: Compute cnt/Eb_tot must use group all

UNSPECIFIED.

E: cnt/Eb_tot is allowed only with cnt pair style

Use cnt pair style.

*/