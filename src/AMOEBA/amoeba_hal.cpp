/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "pair_amoeba.h"
#include <cmath>
#include "atom.h"
#include "neigh_list.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{VDWL,REPULSE,QFER,DISP,MPOLE,POLAR,USOLV,DISP_LONG,MPOLE_LONG,POLAR_LONG};

/* ----------------------------------------------------------------------
   hal = buffered 14-7 Vdwl interactions
   adapted from Tinker ehal1c() routine
------------------------------------------------------------------------- */

void PairAmoeba::hal()
{
  int i,j,ii,jj,itype,jtype,iclass,jclass,iv,jv;
  int special_which;
  double e,de,eps,rdn;
  double fgrp,rv,rv7;
  double xi,yi,zi;
  double xr,yr,zr;
  double redi,rediv;
  double redj,redjv;
  double dedx,dedy,dedz;
  double rho,rho6,rho7;
  double tau,tau7,scal;
  double s1,s2,t1,t2;
  double dt1drho,dt2drho;
  double dtau,gtau;
  double taper,dtaper;
  double rik,rik2,rik3;
  double rik4,rik5;
  double rik6,rik7;
  double vxx,vyy,vzz;
  double vyx,vzx,vzy;
  double factor_hal;
  double vir[6];

  int inum,jnum;
  int *ilist,*jlist,*numneigh,**firstneigh;

  // set cutoffs and taper coeffs

  choose(VDWL);

  // owned atoms

  double **x = atom->x;
  int nlocal = atom->nlocal;

  // neigh list

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // find van der Waals energy and derivatives via neighbor list

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = amtype[i];
    iclass = amtype2class[itype];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    redi = kred[iclass];
    rediv = 1.0 - redi;
    xi = xred[i][0];
    yi = xred[i][1];
    zi = xred[i][2];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      special_which = sbmask15(j);
      factor_hal = special_hal[special_which];
      if (factor_hal == 0.0) continue;
      j &= NEIGHMASK15;
      
      xr = xi - xred[j][0];
      yr = yi - xred[j][1];
      zr = zi - xred[j][2];
      rik2 = xr*xr + yr*yr + zr*zr;

      if (rik2 > off2) continue;

      // compute the energy contribution for this interaction

      jtype = amtype[j];
      jclass = amtype2class[jtype];

      // check for an interaction distance less than the cutoff
      // special_which = 3 is a 1-4 neighbor with its own sigma,epsilon

      rik = sqrt(rik2);
      rv = radmin[iclass][jclass];
      eps = epsilon[iclass][jclass];
      if (special_which == 3) {
        rv = radmin4[iclass][jclass];
        eps = epsilon4[iclass][jclass];
      }
      eps *= factor_hal;

      rv7 = pow(rv,7.0);
      rik6 = pow(rik2,3.0);
      rik7 = rik6 * rik;
      rho = rik7 + ghal*rv7;
      tau = (dhal+1.0) / (rik + dhal*rv);
      tau7 = pow(tau,7.0);
      dtau = tau / (dhal+1.0);
      gtau = eps*tau7*rik6*(ghal+1.0)*pow(rv7/rho,2.0);
      e = eps*tau7*rv7*((ghal+1.0)*rv7/rho-2.0);
      de = -7.0 * (dtau*e+gtau);

      // use energy switching if near the cutoff distance

      if (rik2 > cut2) {
        rik3 = rik2 * rik;
        rik4 = rik2 * rik2;
        rik5 = rik2 * rik3;
        taper = c5*rik5 + c4*rik4 + c3*rik3 + c2*rik2 + c1*rik + c0;
        dtaper = 5.0*c5*rik4 + 4.0*c4*rik3 + 3.0*c3*rik2 + 2.0*c2*rik + c1;
        de = e*dtaper + de*taper;
        e *= taper;
      }

      ehal += e;
      
      // find the chain rule terms for derivative components

      de = de / rik;
      dedx = de * xr;
      dedy = de * yr;
      dedz = de * zr;

      // increment the total van der Waals energy and derivatives
      // if jv < 0, trigger an error, needed H-bond partner is missing
      
      iv = ired2local[i];
      jv = ired2local[j];
      if (jv < 0)
	error->one(FLERR,"AMOEBA hal cannot find H bond partner - "
		   "ghost comm is too short");
      
      if (i == iv) {
        fhal[i][0] += dedx;
        fhal[i][1] += dedy;
        fhal[i][2] += dedz;
      } else {
        fhal[i][0] += dedx*redi;
        fhal[i][1] += dedy*redi;
        fhal[i][2] += dedz*redi;
        fhal[iv][0] += dedx*rediv;
        fhal[iv][1] += dedy*rediv;
        fhal[iv][2] += dedz*rediv;
      }

      if (j == jv) {
        fhal[j][0] -= dedx;
        fhal[j][1] -= dedy;
        fhal[j][2] -= dedz;
      } else {
        redj = kred[jclass];
        redjv = 1.0 - redj;
        fhal[j][0] -= dedx*redj;
        fhal[j][1] -= dedy*redj;
        fhal[j][2] -= dedz*redj;
        fhal[jv][0] -= dedx*redjv;
        fhal[jv][1] -= dedy*redjv;
        fhal[jv][2] -= dedz*redjv;
      }

      // increment the internal virial tensor components

      vxx = xr * dedx;
      vyx = yr * dedx;
      vzx = zr * dedx;
      vyy = yr * dedy;
      vzy = zr * dedy;
      vzz = zr * dedz;

      virhal[0] += vxx;
      virhal[1] += vyy;
      virhal[2] += vzz;
      virhal[3] += vyx;
      virhal[4] += vzx;
      virhal[5] += vzy;

      // energy = e
      // virial = 6-vec vir
      // NOTE: add tally function
      
      if (evflag) {
      }
    }
  }
}
