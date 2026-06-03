// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_wall_srd_region.h"

#include "domain.h"
#include "error.h"
#include "modify.h"
#include "region.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixWallSRDRegion::FixWallSRDRegion(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), region(nullptr), idregion(nullptr), varflag(0), phantom(1),
    force_flag(0)
{
  // syntax: fix ID group-ID wall/srd/region region-ID [phantom yes|no]
  // group-ID is ignored (mirrors fix wall/srd convention -- the SRD group
  // is set by fix srd itself).
  // phantom default = yes (Lamura/Gompper virtual-particle fill);
  // phantom no = bare bounce-back, useful for diagnosing the wall artifact.

  if (narg < 4) error->all(FLERR, "Illegal fix wall/srd/region command");

  idregion = utils::strdup(arg[3]);

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "phantom") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix wall/srd/region command");
      phantom = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else {
      error->all(FLERR, "Unknown keyword for fix wall/srd/region: {}", arg[iarg]);
    }
  }

  // resolve region pointer lazily in init() (region may be defined later
  // in the input script than the fix). We do an early existence check
  // here so syntax errors are caught immediately.

  auto *r = domain->get_region_by_id(idregion);
  if (!r)
    error->all(FLERR, "Region {} for fix wall/srd/region does not exist", idregion);

  // varflag = 1 if region geometry/position can change at runtime
  // (Phase 3: dynamic regions are now allowed; wall_params drives the
  // per-step state update via region->prematch + set_velocity).
  varflag = (r->dynamic || r->varshape) ? 1 : 0;

  // set up output: 3 components of force on this wall, accessible as
  // f_id[1], f_id[2], f_id[3]. Mirrors fix wall/srd's compute_array.

  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extvector = 1;

  fwall[0] = fwall[1] = fwall[2] = 0.0;
  fwall_all[0] = fwall_all[1] = fwall_all[2] = 0.0;
}

/* ---------------------------------------------------------------------- */

FixWallSRDRegion::~FixWallSRDRegion()
{
  delete[] idregion;
}

/* ---------------------------------------------------------------------- */

int FixWallSRDRegion::setmask()
{
  // FixSRD drives everything; we contribute no callbacks of our own.
  return 0;
}

/* ---------------------------------------------------------------------- */

void FixWallSRDRegion::init()
{
  // ensure a matching fix srd exists -- nothing we do is useful otherwise.

  int flag = 0;
  for (int m = 0; m < modify->nfix; m++)
    if (utils::strmatch(modify->fix[m]->style, "^srd")) flag = 1;
  if (!flag) error->all(FLERR, "Cannot use fix wall/srd/region without fix srd");

  // re-resolve region (may have changed across run commands).

  region = domain->get_region_by_id(idregion);
  if (!region)
    error->all(FLERR, "Region {} for fix wall/srd/region does not exist", idregion);
}

/* ----------------------------------------------------------------------
   return component j of net force on this wall (j = 0,1,2 for x,y,z)
------------------------------------------------------------------------- */

double FixWallSRDRegion::compute_vector(int j)
{
  // sum across procs once per access cycle (force_flag reset by wall_params).
  if (!force_flag) {
    MPI_Allreduce(fwall, fwall_all, 3, MPI_DOUBLE, MPI_SUM, world);
    force_flag = 1;
  }
  return fwall_all[j];
}

/* ----------------------------------------------------------------------
   Called by FixSRD.
   flag = 1 at biglist setup time (reneighbor steps).
   flag = 0 at every collision step.
   Phase 1: just zero the per-step force accumulator. In later phases this
   will also drive region->prematch()/set_velocity() for dynamic regions.
------------------------------------------------------------------------- */

void FixWallSRDRegion::wall_params(int /*flag*/)
{
  // For dynamic regions: refresh time-dependent geometry (varshape) and
  // current translational/angular/radial wall velocities. Must run BEFORE
  // any collision check or virtual-particle injection that uses the
  // region's state. fix_wall_gran_region does the same in its post_force.
  if (region->dynamic_check()) {
    region->prematch();
    region->set_velocity();
  }
  fwall[0] = fwall[1] = fwall[2] = 0.0;
  force_flag = 0;
}
