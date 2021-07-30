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

#include "compute_grid.h"
#include "compute_sna_grid.h"
#include "sna.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "force.h"
#include "pair.h"
#include "domain.h"
#include "comm.h"
#include "memory.h"
#include "error.h"
#include "tokenizer.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

ComputeSNAGrid::ComputeSNAGrid(LAMMPS *lmp, int narg, char **arg) :
  ComputeGrid(lmp, narg, arg), cutsq(nullptr),
  radelem(nullptr), wjelem(nullptr)
{

  double rfac0, rmin0;
  int twojmax, switchflag, bzeroflag, bnormflag, wselfallflag;

  // skip over arguments used by base class
  // so that argument positions are identical to 
  // regular per-atom compute

  arg += nargbase;
  narg -= nargbase;

  int ntypes = atom->ntypes;
  int nargmin = 6+2*ntypes;

  if (narg < nargmin) error->all(FLERR,"Illegal compute sna/grid command");

  // default values

  rmin0 = 0.0;
  switchflag = 1;
  bzeroflag = 1;
  quadraticflag = 0;
  chemflag = 0;
  bnormflag = 0;
  wselfallflag = 0;
  nelements = 1;
  
  // process required arguments

  memory->create(radelem,ntypes+1,"sna/grid:radelem"); // offset by 1 to match up with types
  memory->create(wjelem,ntypes+1,"sna/grid:wjelem");

  rcutfac = atof(arg[3]);
  rfac0 = atof(arg[4]);
  twojmax = atoi(arg[5]);

  for(int i = 0; i < ntypes; i++)
    radelem[i+1] = atof(arg[6+i]);
  for(int i = 0; i < ntypes; i++)
    wjelem[i+1] = atof(arg[6+ntypes+i]);

  // construct cutsq

  double cut;
  cutmax = 0.0;
  memory->create(cutsq,ntypes+1,ntypes+1,"sna/grid:cutsq");
  for(int i = 1; i <= ntypes; i++) {
    cut = 2.0*radelem[i]*rcutfac;
    if (cut > cutmax) cutmax = cut;
    cutsq[i][i] = cut*cut;
    for(int j = i+1; j <= ntypes; j++) {
      cut = (radelem[i]+radelem[j])*rcutfac;
      cutsq[i][j] = cutsq[j][i] = cut*cut;
    }
  }

  // process optional args

  int iarg = nargmin;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"rmin0") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute sna/grid command");
      rmin0 = atof(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"switchflag") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute sna/grid command");
      switchflag = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"bzeroflag") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute sna/grid command");
      bzeroflag = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"quadraticflag") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute sna/grid command");
      quadraticflag = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"chem") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute sna/grid command");
      chemflag = 1;
      memory->create(map,ntypes+1,"compute_sna_grid:map");
      nelements = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      for (int i = 0; i < ntypes; i++) {
        int jelem = utils::inumeric(FLERR,arg[iarg+2+i],false,lmp);
        if (jelem < 0 || jelem >= nelements)
          error->all(FLERR,"Illegal compute sna/grid command");
        map[i+1] = jelem;
      }
      iarg += 2+ntypes;
    } else if (strcmp(arg[iarg],"bnormflag") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute sna/grid command");
      bnormflag = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"wselfallflag") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute sna/grid command");
      wselfallflag = atoi(arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal compute sna/grid command");

  }

  snaptr = new SNA(lmp, rfac0, twojmax,
                   rmin0, switchflag, bzeroflag,
                   chemflag, bnormflag, wselfallflag, nelements);

  ncoeff = snaptr->ncoeff;
  nvalues = ncoeff;
  if (quadraticflag) nvalues += (ncoeff*(ncoeff+1))/2;
  size_array_cols = size_array_cols_base + nvalues;
  array_flag = 1;
}

/* ---------------------------------------------------------------------- */

ComputeSNAGrid::~ComputeSNAGrid()
{
  memory->destroy(radelem);
  memory->destroy(wjelem);
  memory->destroy(cutsq);
  delete snaptr;

  if (chemflag) memory->destroy(map);
}

/* ---------------------------------------------------------------------- */

void ComputeSNAGrid::init()
{
  if (force->pair == nullptr)
    error->all(FLERR,"Compute sna/grid requires a pair style be defined");

  if (cutmax > force->pair->cutforce)
    error->all(FLERR,"Compute sna/grid cutoff is longer than pairwise cutoff");

  // need an occasional full neighbor list

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 1;

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"sna/grid") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute sna/grid");
  snaptr->init();
}

/* ---------------------------------------------------------------------- */

void ComputeSNAGrid::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeSNAGrid::compute_array()
{
  invoked_array = update->ntimestep;

  // // invoke full neighbor list (will copy or build if necessary)

  // neighbor->build_one(list);

  // int mynxlo = static_cast<int> (comm->xsplit[comm->myloc[0]] * nx);
  // int mynxhi = static_cast<int> (comm->xsplit[comm->myloc[0]+1] * nx) - 1;

  // int mynylo = static_cast<int> (comm->ysplit[comm->myloc[1]] * ny);
  // int mynyhi = static_cast<int> (comm->ysplit[comm->myloc[1]+1] * ny) - 1;

  // int mynzlo = static_cast<int> (comm->zsplit[comm->myloc[2]] * nz);
  // int mynzhi = static_cast<int> (comm->zsplit[comm->myloc[2]+1] * nz) - 1;

  // int myngridlocal = (nxhi - nxlo + 1) * (nyhi - nylo + 1) * (nzhi - nzlo + 1);

  // compute sna for each gridpoint

  double** const x = atom->x;
  const int* const mask = atom->mask;
  int * const type = atom->type;
  const int ntotal = atom->nlocal + atom->nghost;

  // insure rij, inside, and typej are of size jnum
  
  snaptr->grow_rij(ntotal);

  for (int iz = nzlo; iz <= nzhi; iz++)
    for (int iy = nylo; iy <= nyhi; iy++)
      for (int ix = nxlo; ix <= nxhi; ix++) {
	const int igrid = iz*(nx*ny) + iy*nx + ix;
	const double xtmp = grid[igrid][0];
	const double ytmp = grid[igrid][1];
	const double ztmp = grid[igrid][2];

	// currently, all grid points are type 1
	
	const int itype = 1;
	int ielem = 0;
	if (chemflag)
	  ielem = map[itype];
	const double radi = radelem[itype];

	// rij[][3] = displacements between atom I and those neighbors
	// inside = indices of neighbors of I within cutoff
	// typej = types of neighbors of I within cutoff

	int ninside = 0;
	for (int j = 0; j < ntotal; j++) {

	  // check that j is in compute group

	  if (!(mask[j] & groupbit)) continue;

	  const double delx = xtmp - x[j][0];
	  const double dely = ytmp - x[j][1];
	  const double delz = ztmp - x[j][2];
	  const double rsq = delx*delx + dely*dely + delz*delz;
	  int jtype = type[j];
	  int jelem = 0;
	  if (chemflag)
	    jelem = map[jtype];
	  if (rsq < cutsq[jtype][jtype] && rsq > 1e-20) {
	    snaptr->rij[ninside][0] = delx;
	    snaptr->rij[ninside][1] = dely;
	    snaptr->rij[ninside][2] = delz;
	    snaptr->inside[ninside] = j;
	    snaptr->wj[ninside] = wjelem[jtype];
	    snaptr->rcutij[ninside] = 2.0*radelem[jtype]*rcutfac;
	    snaptr->element[ninside] = jelem; // element index for multi-element snap
	    ninside++;
	  }
	}

	snaptr->compute_ui(ninside, ielem);
	snaptr->compute_zi();
	snaptr->compute_bi(ielem);

      // linear contributions

	for (int icoeff = 0; icoeff < ncoeff; icoeff++)
	  gridlocal[size_array_cols_base+icoeff][iz][iy][ix] = snaptr->blist[icoeff];
	
	// quadratic contributions

	if (quadraticflag) {
	  int ncount = ncoeff;
	  for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
	    double bveci = snaptr->blist[icoeff];
	    gridlocal[size_array_cols_base+ncount++][iz][iy][ix] = 0.5*bveci*bveci;
	    for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++)
	      gridlocal[size_array_cols_base+ncount++][iz][iy][ix] = bveci*snaptr->blist[jcoeff];
	  }
	}
      }
  
  for (int iz = nzlo; iz <= nzhi; iz++)
    for (int iy = nylo; iy <= nyhi; iy++)
      for (int ix = nxlo; ix <= nxhi; ix++) {
	const int igrid = iz*(nx*ny) + iy*nx + ix;
	for (int j = 0; j < nvalues; j++)
	  grid[igrid][size_array_cols_base + j] = gridlocal[size_array_cols_base + j][iz][iy][ix];
      }
  MPI_Allreduce(&grid[0][0],&gridall[0][0],ngrid*size_array_cols,MPI_DOUBLE,MPI_SUM,world);

}


/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double ComputeSNAGrid::memory_usage()
{
  double nbytes = snaptr->memory_usage();    // SNA object
  int n = atom->ntypes+1;
  nbytes += (double)n*sizeof(int);            // map

  return nbytes;
}

