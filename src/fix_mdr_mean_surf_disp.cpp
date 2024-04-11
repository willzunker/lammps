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

#include "fix_mdr_mean_surf_disp.h"

#include "atom.h"
#include "error.h"
#include "force.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "variable.h"
#include "fix_neigh_history.h"
#include "pair.h"
#include "pair_granular.h"
#include "granular_model.h"
#include "neigh_list.h"
#include "region.h"
#include "fix_wall_gran_region.h"
#include <iostream>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace Granular_NS;

/* ---------------------------------------------------------------------- */

FixMDRmeanSurfDisp::FixMDRmeanSurfDisp(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg)
{
 // nothing to initialize
}

// FOR MDR

int FixMDRmeanSurfDisp::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

void FixMDRmeanSurfDisp::setup(int /*vflag*/)
{
  pre_force(0);
}

/* ---------------------------------------------------------------------- */

void FixMDRmeanSurfDisp::pre_force(int)
{
  int tmp1, tmp2;
  int index_Acon0 = atom->find_custom("Acon0",tmp1,tmp2);                 
  int index_ddelta_bar = atom->find_custom("ddelta_bar",tmp1,tmp2);             
  double * Acon0 = atom->dvector[index_Acon0]; 
  double * ddelta_bar = atom->dvector[index_ddelta_bar];

  FixNeighHistory * fix_history = dynamic_cast<FixNeighHistory *>(modify->get_fix_by_id("NEIGH_HISTORY_GRANULAR"));
  PairGranular * pair = dynamic_cast<PairGranular *>(force->pair_match("granular",1));
  NeighList * list = pair->list;
  const int size_history = pair->get_size_history();

  {
  int i,j,k,ii,jj,inum,jnum,itype,jtype;

  int *ilist,*jlist,*numneigh,**firstneigh;
  int *touch,**firsttouch;
  double *history,*allhistory,**firsthistory;

  bool touchflag = false;

  class GranularModel* model;
  class GranularModel** models_list = pair->models_list;
  int ** types_indices = pair->types_indices;

  double **x = atom->x;
  int *type = atom->type;
  double *radius = atom->radius;
  int nlocal = atom->nlocal;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  firsttouch = fix_history->firstflag;
  firsthistory = fix_history->firstvalue;


  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];
    touch = firsttouch[i];
    allhistory = firsthistory[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      jtype = type[j];
      model = models_list[types_indices[itype][jtype]];

      // Reset model and copy initial geometric data
      model->xi = x[i];
      model->xj = x[j];
      model->radi = radius[i];
      model->radj = radius[j];
      model->i = i;
      model->j = j;

      model->touch = touch[jj];

      touchflag = model->check_contact();

      if (!touchflag) {
        touch[jj] = 0;
        history = &allhistory[size_history * jj];
        for (k = 0; k < size_history; k++) history[k] = 0.0;
        continue;
      }

      touch[jj] = 1;


      history = &allhistory[size_history * jj];
      model->history = history;


      const double wij = 1.0;
      const double delta = model->radsum - sqrt(model->rsq);

      if (Acon0[j] != 0.0) {
        const double delta_offset0 = history[0];
        const double ddelta = delta/2.0 - delta_offset0; // Divide by 2.0 since we are storing 1/2 deltan in main MDR script
        const double Ac_offset0 = history[18];
        ddelta_bar[j] += wij*Ac_offset0/Acon0[j]*ddelta; // Multiply by 0.5 since displacement is shared equally between deformable particles.
      }

      if (Acon0[i] != 0.0) {
        const double delta_offset1 = history[1];
        const double ddelta = delta/2.0 - delta_offset1; // Divide by 2.0 since we are storing 1/2 deltan in main MDR script
        const double Ac_offset1 = history[19];
        ddelta_bar[i] += wij*Ac_offset1/Acon0[i]*ddelta;
      }

    }
  }
}

  auto fix_list = modify->get_fix_by_style("wall/gran/region");

  for (int w = 0; w < fix_list.size(); w++) {

    FixWallGranRegion* fix = dynamic_cast<FixWallGranRegion*>(fix_list[w]);
    GranularModel * model = fix->model;
    Region * region = fix->region;

    {
    int i, m, nc, iwall;
    double vwall[3];
    bool touchflag = false;

    int history_update = 1;
    model->history_update = history_update;

    int regiondynamic = region->dynamic_check();
    if (!regiondynamic) vwall[0] = vwall[1] = vwall[2] = 0.0;

    double **x = atom->x;
    double *radius = atom->radius;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    if (regiondynamic) {
      region->prematch();
      region->set_velocity();
    }

    if (fix->peratom_flag) fix->clear_stored_contacts();

    model->radj = 0.0;

    for (i = 0; i < nlocal; i++) {
      if (!(mask[i] & groupbit)) continue;
      if (! region->match(x[i][0], x[i][1], x[i][2])) continue;

      nc = region->surface(x[i][0], x[i][1], x[i][2], radius[i] + model->pulloff_distance(radius[i], 0.0));

        if (nc == 0) {
          fix->ncontact[i] = 0;
          continue;
        }
        if (nc == 1) {
          fix->c2r[0] = 0;
          iwall = region->contact[0].iwall;
          if (fix->ncontact[i] == 0) {
            fix->ncontact[i] = 1;
            fix->walls[i][0] = iwall;
            for (m = 0; m < size_history; m++) fix->history_many[i][0][m] = 0.0;
          } else if (fix->ncontact[i] > 1 || iwall != fix->walls[i][0])
            fix->update_contacts(i, nc);
        } else
          fix->update_contacts(i, nc);


      // process current contacts
      for (int ic = 0; ic < nc; ic++) {

        // Reset model and copy initial geometric data
        model->dx[0] = region->contact[ic].delx;
        model->dx[1] = region->contact[ic].dely;
        model->dx[2] = region->contact[ic].delz;
        model->radi = radius[i];
        model->radj = region->contact[ic].radius;
        model->r = region->contact[ic].r;

        if (model->beyond_contact) model->touch = fix->history_many[i][fix->c2r[ic]][0];

        touchflag = model->check_contact();

        const double wij = 1.0;

        if (Acon0[i] != 0.0) {
          const double delta = model->radsum - model->r;
          const double delta_offset0 = fix->history_many[i][fix->c2r[ic]][0];
          const double ddelta = delta - delta_offset0; 
          const double Ac_offset0 = fix->history_many[i][fix->c2r[ic]][18];
          ddelta_bar[i] += wij*Ac_offset0/Acon0[i]*ddelta; // Multiply by 0.5 since displacement is shared equally between deformable particles.
          //std::cout << delta << ", " << delta_offset0 << " || " << Ac_offset0 << ", " << Acon0[i] << ", " << ddelta << ", " << ddelta_bar[i] << std::endl;
        }
      }
    } 
    }
  }

}

//std::cout << radius[i] << ", " << dR << ", " << dRnumerator[i] << ", " << dRdenominator[i] << ", " << dRdenominator[i] - 4.0*M_PI*pow(R,2.0)  << std::endl;
//std::cout << "Fix radius update setup has been entered !!!" << std::endl;
//std::cout << Ro[0] << ", " << Vgeo[0] << ", " << Velas[0] << std::endl;