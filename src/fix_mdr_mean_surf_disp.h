/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(mdr/mean/surf/disp,FixMDRmeanSurfDisp);
// clang-format on
#else

#ifndef LMP_FIX_MDR_MEAN_SURF_DISP_H
#define LMP_FIX_MDR_MEAN_SURF_DISP_H

#include "fix.h"

namespace LAMMPS_NS {

class FixMDRmeanSurfDisp : public Fix {
 public:
  FixMDRmeanSurfDisp(class LAMMPS *, int, char **);
  int setmask() override;
  void setup(int) override;
  void pre_force(int) override; // FOR MDR

 private:

};

}    // namespace LAMMPS_NS

#endif
#endif
