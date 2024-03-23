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

#include "fix_mdr_radius_update.h"

#include "atom.h"
#include "error.h"
#include "force.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "variable.h"
#include <iostream>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixMDRradiusUpdate::FixMDRradiusUpdate(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg)
{
 // nothing to initialize
}

// FOR MDR

int FixMDRradiusUpdate::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

void FixMDRradiusUpdate::setup(int /*vflag*/)
{
  
  int tmp1, tmp2;
  int index_Ro = atom->find_custom("Ro",tmp1,tmp2);
  int index_Rold = atom->find_custom("Rold",tmp1,tmp2);
  double * Ro = atom->dvector[index_Ro];
  double * Rold = atom->dvector[index_Rold];

  double *radius = atom->radius;
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    Ro[i] = radius[i];
    Rold[i] = radius[i];
  }

  std::cout << "Fix radius update setup has been entered !!!" << std::endl;
  std::cout << Ro[0] << ", " << Rold[0] << std::endl;

  end_of_step();
}

/* ---------------------------------------------------------------------- */

// FOR MDR, DO WHATEVER YOUR FIX NEEDS TO DO.

void FixMDRradiusUpdate::end_of_step()
{
  double dR = 0.1;
  double *radius = atom->radius;
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    radius[i] += dR;
  }
}

