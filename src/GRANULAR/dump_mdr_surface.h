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

#ifdef DUMP_CLASS
// clang-format off
DumpStyle(mdr/surface,DumpMDRSurface);
// clang-format on
#else

#ifndef LMP_DUMP_MDR_SURFACE_H
#define LMP_DUMP_MDR_SURFACE_H

#include "dump.h"

#include <map>
#include <string>
#include <vector>

namespace LAMMPS_NS {

namespace Granular_NS {
  class GranSubModNormalMDR;
}

class DumpMDRSurface : public Dump {
 public:
  DumpMDRSurface(class LAMMPS *, int, char **);
  ~DumpMDRSurface() override;
  void write() override;

 protected:
  void init_style() override;
  void write_header(bigint) override {}
  void pack(tagint *) override {}
  void write_data(int, double *) override {}

  int ntheta, nphi;          // surface mesh resolution
  int csv_flag;              // also write raw per-contact CSV (debug / MATLAB compat)
  char *csvfile;             // template for csv output (with '*')

  // resolved on init / first write
  class PairGranular *pair;
  Granular_NS::GranSubModNormalMDR *mdr_model;
  class FixNeighHistory *fix_hist;
  int index_Ro, index_psi, index_sigmaxx, index_sigmayy, index_sigmazz;

  // persistent plastic imprints, keyed by atom tag: the deepest cap params seen
  // for each plastic contact, so a contact that later separates still shows its
  // residual (elastically sprung-back) flat instead of recovering to a sphere
  struct PersistCap {
    double n[3];
    double R, delta, deltae, A, B, a_na, deltamax, deltaR, amax;
  };
  std::map<tagint, std::vector<PersistCap>> persist;

  // pvd time-series bookkeeping (proc 0 only)
  std::string pvdname;
  std::vector<std::pair<double, std::string>> pvd_entries;

  void find_mdr_model();
  std::string this_step_name(const char *templatename);
  void write_pvd();
};

}    // namespace LAMMPS_NS

#endif
#endif
