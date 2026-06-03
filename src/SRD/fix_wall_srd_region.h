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
FixStyle(wall/srd/region,FixWallSRDRegion);
// clang-format on
#else

#ifndef LMP_FIX_WALL_SRD_REGION_H
#define LMP_FIX_WALL_SRD_REGION_H

#include "fix.h"

namespace LAMMPS_NS {

// Phase 1: SRD wall whose geometry is defined by a LAMMPS region.
// Mirrors the role of FixWallSRD for the planar-face case, but delegates
// all geometry (point-in-region, surface distance, surface normal,
// surface velocity) to the Region object. Allows cylinder, sphere,
// intersect, union, etc. -- anything region supports.
//
// Phase 1 limitations:
//   - Static regions only (no motion / no time-varying shape).
//   - Inexact collision only (push-to-surface at end of step).
//   - Force accumulator is a single 3-vector per fix (no per-iwall split
//     for compound regions yet).

class FixWallSRDRegion : public Fix {
 public:
  class Region *region;      // pointer to the region (resolved in init)
  char *idregion;            // region id string

  int varflag;               // 1 if region is dynamic / time-varying
  int phantom;               // 1 = Lamura/Gompper virtual-particle fill (default)
                             // 0 = bare bounce-back (diagnostic / comparison)
  double fwall[3];           // momentum/force accumulated on this wall

  FixWallSRDRegion(class LAMMPS *, int, char **);
  ~FixWallSRDRegion() override;
  int setmask() override;
  void init() override;
  double compute_vector(int) override;

  // Called by FixSRD: flag=1 at biglist setup (reneighbor), flag=0 at
  // every collision step. Zeros forces; in later phases will also update
  // region motion state via region->prematch()/set_velocity().
  void wall_params(int flag);

 private:
  double fwall_all[3];
  int force_flag;
};

}    // namespace LAMMPS_NS

#endif
#endif
