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
  //std::cout << "New Step" << std::endl;

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
  int i,j,k,lv1,ii,jj,inum,jnum,itype,jtype,ktype;

  int *ilist,*jlist,*numneigh,**firstneigh;
  int *touch,**firsttouch;
  double *history_ij,*history_ik,*history_jk,*history_kj,*allhistory,*allhistory_j,*allhistory_k,**firsthistory;

  bool touchflag = false;

  //class GranularModel* model;
  //class GranularModel** models_list = pair->models_list;
  //int ** types_indices = pair->types_indices;

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

  // contact penalty calculation
  for (int ii = 0; ii < inum; ii++) {
    const int i = ilist[ii];
    const double xtmp = x[i][0];
    const double ytmp = x[i][1];
    const double ztmp = x[i][2];
    allhistory = firsthistory[i];
    double radi = radius[i]; 
    jlist = firstneigh[i];
    jnum = numneigh[i];
      for (int jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;
        jtype = type[j];
        double radj = radius[j];        
        const double delx_ij = x[j][0] - xtmp;
        const double dely_ij = x[j][1] - ytmp;
        const double delz_ij = x[j][2] - ztmp;
        const double rsq_ij = delx_ij*delx_ij + dely_ij*dely_ij + delz_ij*delz_ij;
        const double r_ij = sqrt(rsq_ij);
        const double rinv_ij = 1.0/r_ij;
        const double radsum_ij = radi + radj;
        const double deltan_ij = radsum_ij - r_ij;
        if (deltan_ij >= 0.0) {
          for (int kk = jj; kk < jnum; kk++) {
            k = jlist[kk]; 
            k &= NEIGHMASK;
            ktype = type[k];
            if (kk != jj) {
              const double delx_ik = x[k][0] - xtmp;
              const double dely_ik = x[k][1] - ytmp;
              const double delz_ik = x[k][2] - ztmp;
              const double rsq_ik = delx_ik*delx_ik + dely_ik*dely_ik + delz_ik*delz_ik;
              const double r_ik = sqrt(rsq_ik);
              const double rinv_ik = 1.0/r_ik;
              const double radk = radius[k];
              const double radsum_ik = radi + radk;
              const double deltan_ik = radsum_ik - r_ik;
              const double delx_jk = x[k][0] - x[j][0];
              const double dely_jk = x[k][1] - x[j][1];
              const double delz_jk = x[k][2] - x[j][2];
              const double rsq_jk = delx_jk*delx_jk + dely_jk*dely_jk + delz_jk*delz_jk;
              const double r_jk = sqrt(rsq_jk);
              const double rinv_jk = 1.0/r_jk;
              const double radsum_jk = radj + radk;
              const double deltan_jk = radsum_jk - r_jk;
              if (deltan_ik >= 0.0 && deltan_jk >= 0.0) {
                
                // pull ij history
                history_ij = &allhistory[size_history * jj];
                double * pij = &history_ij[22]; // penalty for contact i and j

                // pull ik history
                history_ik = &allhistory[size_history * kk];
                double * pik = &history_ik[22]; // penalty for contact i and k

                // we don't know if who owns the contact ahead of time, k might be in j's neigbor list or vice versa, so we need to manually search to figure out the owner
                // check if k is in the neighbor list of j
                double * pjk = NULL; 
                int * const jklist = firstneigh[j];
                const int jknum = numneigh[j];
                for (int jk = 0; jk < jknum; jk++) {
                  const int kneigh = jklist[jk] & NEIGHMASK;
                  if (k == kneigh) {
                    allhistory_j = firsthistory[j];
                    history_jk = &allhistory_j[size_history * jk];
                    pjk = &history_jk[22]; // penalty for contact j and k
                    break;
                  }
                }

                // check if j is in the neighbor list of k
                if (pjk == NULL) {
                  int * const kjlist = firstneigh[k];
                  const int kjnum = numneigh[k];
                  for (int kj = 0; kj < kjnum; kj++) {
                    const int jneigh = kjlist[kj] & NEIGHMASK;
                    if (j == jneigh) {
                      allhistory_k = firsthistory[k];
                      history_kj = &allhistory_k[size_history * kj];
                      pjk = &history_kj[22]; // penalty for contact j and k
                      break;
                    }
                  }         
                }

                std::vector<double> distances = {r_ij,r_ik,r_jk};
                auto maxElement = std::max_element(distances.begin(), distances.end());
                double maxValue = *maxElement;
                int maxIndex = std::distance(distances.begin(), maxElement);
                if (maxIndex == 0) { // the central particle is k
                  const double enx_ki = -delx_ik * rinv_ik;
                  const double eny_ki = -dely_ik * rinv_ik;
                  const double enz_ki = -delz_ik * rinv_ik;
                  const double enx_kj = -delx_jk * rinv_jk;
                  const double eny_kj = -dely_jk * rinv_jk;
                  const double enz_kj = -delz_jk * rinv_jk;
                  const double alpha = std::acos(enx_ki*enx_kj + eny_ki*eny_kj + enz_ki*enz_kj); 
                  pij[0] += 1.0/( 1.0 + std::exp(-50.0*(alpha/M_PI - 1.0/2.0)) ); 
                } else if (maxIndex == 1) { // the central particle is j
                  const double enx_ji = -delx_ij * rinv_ij;
                  const double eny_ji = -dely_ij * rinv_ij;
                  const double enz_ji = -delz_ij * rinv_ij;
                  const double enx_jk = delx_jk * rinv_jk;
                  const double eny_jk = dely_jk * rinv_jk;
                  const double enz_jk = delz_jk * rinv_jk;
                  const double alpha = std::acos(enx_ji*enx_jk + eny_ji*eny_jk + enz_ji*enz_jk); 
                  pik[0] += 1.0/( 1.0 + std::exp(-50.0*(alpha/M_PI - 1.0/2.0)) );
                } else { // the central particle is i
                  if (j < atom->nlocal || k < atom->nlocal) {
                    const double enx_ij = delx_ij * rinv_ij;
                    const double eny_ij = dely_ij * rinv_ij;
                    const double enz_ij = delz_ij * rinv_ij;
                    const double enx_ik = delx_ik * rinv_ik;
                    const double eny_ik = dely_ik * rinv_ik;
                    const double enz_ik = delz_ik * rinv_ik;
                    const double alpha = std::acos(enx_ij*enx_ik + eny_ij*eny_ik + enz_ij*enz_ik); 
                    pjk[0] += 1.0/( 1.0 + std::exp(-50.0*(alpha/M_PI - 1.0/2.0)) );
                  }
                }
              }
            }
          }
        }
      }
  } 
  }


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

      // is it necessary to clear the history here???
      if (!touchflag) {
        touch[jj] = 0;
        history = &allhistory[size_history * jj];
        for (k = 0; k < size_history; k++) history[k] = 0.0;
        continue;
      }

      touch[jj] = 1;

      history = &allhistory[size_history * jj];
      model->history = history;

      const double delta = model->radsum - sqrt(model->rsq);

      if (Acon0[j] != 0.0) {
        const double delta_offset0 = history[0];
        const double ddelta = delta/2.0 - delta_offset0; // Divide by 2.0 since we are storing 1/2 deltan in main MDR script
        const double Ac_offset0 = history[18];
        ddelta_bar[j] += Ac_offset0/Acon0[j]*ddelta; // Multiply by 0.5 since displacement is shared equally between deformable particles.
      }

      if (Acon0[i] != 0.0) {
        const double delta_offset1 = history[1];
        const double ddelta = delta/2.0 - delta_offset1; // Divide by 2.0 since we are storing 1/2 deltan in main MDR script
        const double Ac_offset1 = history[19];
        ddelta_bar[i] += Ac_offset1/Acon0[i]*ddelta;
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
        }
      }
    } 
    }
  }

}