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
  // assign correct value to initially non-zero MDR particle history variables 
  int tmp1, tmp2;
  int index_Ro = atom->find_custom("Ro",tmp1,tmp2);
  int index_Vgeo = atom->find_custom("Vgeo",tmp1,tmp2);
  int index_Velas = atom->find_custom("Velas",tmp1,tmp2);
  double * Ro = atom->dvector[index_Ro];
  double * Vgeo = atom->dvector[index_Vgeo];
  double * Velas = atom->dvector[index_Velas];

  double *radius = atom->radius;
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    Ro[i] = radius[i];
    Vgeo[i] = 4.0/3.0*M_PI*pow(Ro[i],3.0);
    Velas[i] = 4.0/3.0*M_PI*pow(Ro[i],3.0);
  }
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixMDRradiusUpdate::end_of_step()
{
  // update the apparent radius of every particle

  int tmp1, tmp2;
  int index_Ro = atom->find_custom("Ro",tmp1,tmp2);
  int index_Vgeo = atom->find_custom("Vgeo",tmp1,tmp2);
  int index_Velas = atom->find_custom("Velas",tmp1,tmp2);
  int index_Vcaps = atom->find_custom("Vcaps",tmp1,tmp2);
  int index_eps_bar = atom->find_custom("eps_bar",tmp1,tmp2);
  int index_dRnumerator = atom->find_custom("dRnumerator",tmp1,tmp2);
  int index_dRdenominator = atom->find_custom("dRdenominator",tmp1,tmp2);
  double * Ro = atom->dvector[index_Ro];
  double * Vgeo = atom->dvector[index_Vgeo];
  double * Velas = atom->dvector[index_Velas];
  double * Vcaps = atom->dvector[index_Vcaps];
  double * eps_bar = atom->dvector[index_eps_bar];
  double * dRnumerator = atom->dvector[index_dRnumerator];
  double * dRdenominator = atom->dvector[index_dRdenominator];

  double *radius = atom->radius;
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    const double R = radius[i];
    const double Vo = 4.0/3.0*M_PI*pow(Ro[i],3.0);
    (Vgeo[i] < Vo) ? Vgeo[i] = 4.0/3.0*M_PI*pow(R,3.0) - Vcaps[i] : Vgeo[i] = Vo;

    const double dR = std::max(dRnumerator[i]/(dRdenominator[i] - 4.0*M_PI*pow(R,2.0)),0.0);
    const double psi = 1.0;
    const double psi_b = 0.08;
    if (psi_b < psi) { 
      radius[i] += dR;
      
    }

    Velas[i] = Vo*(1.0 + eps_bar[i]);
    Vcaps[i] = 0.0;
    eps_bar[i] = 0.0;
    dRnumerator[i] = 0.0;
    dRdenominator[i] = 0.0;
  }
}

//std::cout << radius[i] << ", " << dR << ", " << dRnumerator[i] << ", " << dRdenominator[i] << ", " << dRdenominator[i] - 4.0*M_PI*pow(R,2.0)  << std::endl;
//std::cout << "Fix radius update setup has been entered !!!" << std::endl;
//std::cout << Ro[0] << ", " << Vgeo[0] << ", " << Velas[0] << std::endl;