/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
es   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include <string.h>
#include "npair_half_size_multi_newtoff.h"
#include "neigh_list.h"
#include "atom.h"
#include "atom_vec.h"
#include "my_page.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NPairHalfSizeMultiNewtoff::NPairHalfSizeMultiNewtoff(LAMMPS *lmp) : NPair(lmp) {}

/* ----------------------------------------------------------------------
   size particles
   binned neighbor list construction with partial Newton's 3rd law
   multi stencil is igroup-jgroup dependent      
   each owned atom i checks own bin and other bins in stencil
   pair stored once if i,j are both owned and i < j
   pair stored by me if j is ghost (also stored by proc owning j)
------------------------------------------------------------------------- */

void NPairHalfSizeMultiNewtoff::build(NeighList *list)
{
  int i,j,k,n,itype,jtype,igroup,jgroup,ibin,jbin,ns;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  double radi,radsum,cutdistsq;
  int *neighptr,*s;
  int js;

  double **x = atom->x;
  double *radius = atom->radius;
  int *type = atom->type;
  int *mask = atom->mask;
  tagint *molecule = atom->molecule;
  int nlocal = atom->nlocal;
  if (includegroup) nlocal = atom->nfirst;

  int history = list->history;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  MyPage<int> *ipage = list->ipage;

  int mask_history = 3 << SBBITS;

  int inum = 0;
  ipage->reset();

  for (i = 0; i < nlocal; i++) {
    n = 0;
    neighptr = ipage->vget();
    itype = type[i];
    igroup = map_type_multi[itype];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];

    ibin = atom2bin[i];
    
    // loop through stencils for all groups    
    for (jgroup = 0; jgroup < n_multi_groups; jgroup++) {
        
      // if same group use own bin
      if(igroup == jgroup) jbin = ibin;
	  else jbin = coord2bin(x[i], jgroup);
      
      // loop over all atoms in other bins in stencil including self
      // only store pair if i < j
      // stores own/own pairs only once
      // stores own/ghost pairs on both procs      
      // use full stencil for all group combinations

      s = stencil_multi[igroup][jgroup];
      ns = nstencil_multi[igroup][jgroup];
      
      for (k = 0; k < ns; k++) {
	    js = binhead_multi[jgroup][jbin + s[k]];
	    for (j = js; j >=0; j = bins[j]) {
	      if (j <= i) continue;
           
          jtype = type[j];
          if (exclude && exclusion(i,j,itype,jtype,mask,molecule)) continue;
          
	      delx = xtmp - x[j][0];
	      dely = ytmp - x[j][1];
	      delz = ztmp - x[j][2];
	      rsq = delx*delx + dely*dely + delz*delz;
	      radsum = radi + radius[j];
	      cutdistsq = (radsum+skin) * (radsum+skin);

	      if (rsq <= cutdistsq) {
	        if (history && rsq < radsum*radsum) 
	          neighptr[n++] = j ^ mask_history;
	        else
	          neighptr[n++] = j;
	      }
	    }
      }
    }

    ilist[inum++] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    ipage->vgot(n);
    if (ipage->status())
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
  }

  list->inum = inum;
}
