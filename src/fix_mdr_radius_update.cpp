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

/* ----------------------------------------------------------------------
   Contributing authors:
   William Zunker (MIT), Sachith Dunatunga (MIT),
   Dan Bolintineanu (SNL), Joel Clemmer (SNL)
----------------------------------------------------------------------- */

#include "fix_mdr_radius_update.h"

#include "atom.h"
#include "error.h"
#include "force.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "variable.h"
#include <iostream>
#include "csv_writer.h"
#include "granular_model.h"
#include "pair_granular.h"
#include "pair.h"
#include "gran_sub_mod_normal.h"
#include <iomanip> 
#include <sstream>

using namespace LAMMPS_NS;
using namespace Granular_NS;
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
  mask |= PRE_FORCE | END_OF_STEP;
  return mask;
}

void FixMDRradiusUpdate::pre_force(int)
{

  PairGranular * pair = dynamic_cast<PairGranular *>(force->pair_match("granular",1));
  class GranularModel* model;
  class GranularModel** models_list = pair->models_list;
  class GranSubModNormalMDR* norm_model = nullptr;
  for (int i = 0; i < pair->nmodels; i++) {
    model = models_list[i];
    if (model->normal_model->name == "mdr") norm_model = dynamic_cast<GranSubModNormalMDR *>(model->normal_model);
  }
  if (norm_model == nullptr) error->all(FLERR, "Did not find mdr model");

  // assign correct value to initially non-zero MDR particle history variables 
  int tmp1, tmp2;
  int index_Ro = atom->find_custom("Ro",tmp1,tmp2);
  int index_Vgeo = atom->find_custom("Vgeo",tmp1,tmp2);
  int index_Velas = atom->find_custom("Velas",tmp1,tmp2);
  int index_Atot = atom->find_custom("Atot",tmp1,tmp2);                   
  int index_psi = atom->find_custom("psi",tmp1,tmp2);
  int index_psi_b = atom->find_custom("psi_b",tmp1,tmp2);
  int index_sigmaxx = atom->find_custom("sigmaxx",tmp1,tmp2);             
  int index_sigmayy = atom->find_custom("sigmayy",tmp1,tmp2);               
  int index_sigmazz = atom->find_custom("sigmazz",tmp1,tmp2);   
  int index_history_setup_flag = atom->find_custom("history_setup_flag",tmp1,tmp2); 
  int index_contacts = atom->find_custom("contacts",tmp1,tmp2);    
  int index_adhesive_length = atom->find_custom("adhesive_length",tmp1,tmp2);         
  double * Ro = atom->dvector[index_Ro];
  double * Vgeo = atom->dvector[index_Vgeo];
  double * Velas = atom->dvector[index_Velas];
  double * Atot = atom->dvector[index_Atot];
  double * psi = atom->dvector[index_psi];
  double * psi_b = atom->dvector[index_psi_b];
  double * sigmaxx = atom->dvector[index_sigmaxx];
  double * sigmayy = atom->dvector[index_sigmayy];
  double * sigmazz = atom->dvector[index_sigmazz];
  double * history_setup_flag = atom->dvector[index_history_setup_flag];
  double * contacts = atom->dvector[index_contacts];
  double * adhesive_length = atom->dvector[index_adhesive_length];
  
  double *radius = atom->radius;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) { 
    if (history_setup_flag[i] < 1e-16) {
      Ro[i] = radius[i];
      Vgeo[i] = 4.0/3.0*M_PI*pow(Ro[i],3.0);
      Velas[i] = 4.0/3.0*M_PI*pow(Ro[i],3.0);
      Atot[i] = 4.0*M_PI*pow(Ro[i],2.0);
      psi[i] = 1.0;
      psi_b[i] = norm_model->psi_b;
      history_setup_flag[i] = 1.0;
    }
    sigmaxx[i] = 0.0;
    sigmayy[i] = 0.0;
    sigmazz[i] = 0.0;
    contacts[i] = 0.0;
    adhesive_length[i] = 0.0;
  }
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
  int index_Acon0 = atom->find_custom("Acon0",tmp1,tmp2);                 
  int index_Acon1 = atom->find_custom("Acon1",tmp1,tmp2);        
  int index_Atot = atom->find_custom("Atot",tmp1,tmp2);                   
  int index_Atot_sum = atom->find_custom("Atot_sum",tmp1,tmp2);   
  int index_ddelta_bar = atom->find_custom("ddelta_bar",tmp1,tmp2);       
  int index_psi = atom->find_custom("psi",tmp1,tmp2);
  int index_psi_b = atom->find_custom("psi_b",tmp1,tmp2);    
  double * Ro = atom->dvector[index_Ro];
  double * Vgeo = atom->dvector[index_Vgeo];
  double * Velas = atom->dvector[index_Velas];
  double * Vcaps = atom->dvector[index_Vcaps];
  double * eps_bar = atom->dvector[index_eps_bar];
  double * dRnumerator = atom->dvector[index_dRnumerator];
  double * dRdenominator = atom->dvector[index_dRdenominator];
  double * Acon0 = atom->dvector[index_Acon0]; 
  double * Acon1 = atom->dvector[index_Acon1];
  double * Atot = atom->dvector[index_Atot]; 
  double * Atot_sum = atom->dvector[index_Atot_sum];
  double * ddelta_bar = atom->dvector[index_ddelta_bar];
  double * psi = atom->dvector[index_psi];
  double * psi_b = atom->dvector[index_psi_b];

  double *radius = atom->radius;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
  
    const double R = radius[i];
    Atot[i] = 4.0*M_PI*pow(R,2.0) + Atot_sum[i];

    const double Vo = 4.0/3.0*M_PI*pow(Ro[i],3.0);
    const double Vgeoi = 4.0/3.0*M_PI*pow(R,3.0) - Vcaps[i];
    Vgeo[i] = std::min(Vgeoi,Vo);

    const double Afree = Atot[i] - Acon1[i];
    psi[i] = Afree/Atot[i];

    const double dR = std::max(dRnumerator[i]/(dRdenominator[i] - 4.0*M_PI*pow(R,2.0)),0.0);
    if (psi_b[i] < psi[i]) { 
      radius[i] += dR;
    }

    Velas[i] = Vo*(1.0 + eps_bar[i]);
    Vcaps[i] = 0.0;
    eps_bar[i] = 0.0;
    dRnumerator[i] = 0.0;
    dRdenominator[i] = 0.0;
    Acon0[i] = Acon1[i];
    Acon1[i] = 0.0;
    Atot_sum[i] = 0.0;
    ddelta_bar[i] = 0.0;
  }
}