// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_store.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "memory.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

enum{UNKNOWN,GLOBAL,PERATOM};

/* ---------------------------------------------------------------------- */

FixStore::FixStore(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg),
vstore(nullptr), astore(nullptr), rbuf(nullptr)
{
  if (narg != 6 && narg != 7) error->all(FLERR,"Illegal fix store command");

  // 4th arg determines GLOBAL vs PERATOM values
  // global syntax: id group style global n1 n2
  //   N2 = 1 is vector, N2 > 1 is array, no tensor allowed (yet)
  // peratom syntax: id group style peratom 0/1 n1 n2 (n3), last arg optional
  //   0/1 flag = not-store or store peratom values in restart file
  //   N2 = 1 and no n3 is vector, N2 > 1 and no n3 is array, N3 = tensor
  //   nvalues = # of peratom values, N = 1 is vector, N > 1 is array

  disable = 0;
  flavor = UNKNOWN;

  if (strcmp(arg[3],"global") == 0) flavor = GLOBAL;
  else if (strcmp(arg[3],"peratom") == 0) flavor = PERATOM;
  else error->all(FLERR,"Illegal fix store command");

  // GLOBAL values are always written to restart file
  // PERATOM restart_peratom is set by caller

  vecflag = arrayflag = tensorflag = 0;

  if (flavor == GLOBAL) {
    restart_global = 1;
    n1 = utils::inumeric(FLERR,arg[4],false,lmp);
    n2 = utils::inumeric(FLERR,arg[5],false,lmp);
    if (narg == 7) error->all(FLERR,"Illegal fix store command");
    if (n1 <= 0 || n2 <= 0) error->all(FLERR,"Illegal fix store command");
    if (n2 == 1) vecflag = 1;
    else arrayflag = 1;
    nrow = n1;
    ncol = n2;
  }

  if (flavor == PERATOM) {
    restart_peratom = utils::inumeric(FLERR,arg[4],false,lmp);
    n2 = utils::inumeric(FLERR,arg[5],false,lmp);
    if (narg == 7) n3 = utils::inumeric(FLERR,arg[6],false,lmp);
    else n3 = 1;
    if (restart_peratom < 0 || restart_peratom > 1)
      error->all(FLERR,"Illegal fix store command");
    if (n2 <= 0 || n3 <= 0) error->all(FLERR,"Illegal fix store command");
    if (n2 == 1 && narg == 6) vecflag = 1;
    else if (narg == 6) arrayflag = 1;
    else tensorflag = 1;
    nvalues = n2*n3;
    nbytes = n2*n3 * sizeof(double);
  }

  vstore = NULL;
  astore = NULL;
  tstore = NULL;

  // allocate data struct and restart buffer rbuf
  // for PERATOM, register with Atom class

  if (flavor == GLOBAL) {
    if (vecflag) memory->create(vstore,n1,"fix/store:vstore");
    else if (arrayflag) memory->create(astore,n1,n2,"fix/store:astore");
    memory->create(rbuf,n1*n2+2,"fix/store:rbuf");
  }

  if (flavor == PERATOM) {
    grow_arrays(atom->nmax);
    atom->add_callback(Atom::GROW);
    if (restart_peratom) atom->add_callback(Atom::RESTART);
    rbuf = nullptr;
  }

  // zero the storage
  // PERATOM may be comm->exchanged before filled by caller

  if (flavor == GLOBAL) {
    if (vecflag) {
      for (int i = 0; i < n1; i++) 
        vstore[i] = 0.0;
    } else if (arrayflag) {
      for (int i = 0; i < n1; i++)
        for (int j = 0; j < n2; j++)
          astore[i][j] = 0.0;
    }
  }

  if (flavor == PERATOM) {
    int nlocal = atom->nlocal;
    if (vecflag) {
      for (int i = 0; i < nlocal; i++) 
        vstore[i] = 0.0;
    } else if (arrayflag) {
      for (int i = 0; i < nlocal; i++)
        for (int j = 0; j < n2; j++)
          astore[i][j] = 0.0;
    } else if (tensorflag) {
      for (int i = 0; i < nlocal; i++)
        for (int j = 0; j < n2; j++)
          for (int k = 0; k < n3; k++)
            tstore[i][j][k] = 0.0;
    }
    maxexchange = nvalues;
  }
}

/* ---------------------------------------------------------------------- */

FixStore::~FixStore()
{
  // unregister callbacks to this fix from Atom class

  if (flavor == PERATOM) {
    atom->delete_callback(id,Atom::GROW);
    if (restart_peratom) atom->delete_callback(id,Atom::RESTART);
  }

  memory->destroy(vstore);
  memory->destroy(astore);
  memory->destroy(tstore);
  memory->destroy(rbuf);
}

/* ---------------------------------------------------------------------- */

int FixStore::setmask()
{
  int mask = 0;
  return mask;
}

/* ----------------------------------------------------------------------
   reset size of global vector/array
   invoked by caller if size is unknown at time this fix is instantiated
   caller will do subsequent initialization
------------------------------------------------------------------------- */

 void FixStore::reset_global(int n1_caller, int n2_caller)
{
  memory->destroy(vstore);
  memory->destroy(astore);
  memory->destroy(rbuf);
  vstore = nullptr;
  astore = nullptr;

  vecflag = arrayflag = 0;
  if (n2_caller == 1) vecflag = 1;
  else arrayflag = 1;

  n1 = n1_caller;
  n2 = n2_caller;
  if (vecflag) memory->create(vstore,n1,"fix/store:vstore");
  else if (arrayflag) memory->create(astore,n1,n2,"fix/store:astore");
  memory->create(rbuf,n1*n2+2,"fix/store:rbuf");
}

/* ----------------------------------------------------------------------
   write global vector/array to restart file
------------------------------------------------------------------------- */

void FixStore::write_restart(FILE *fp)
{
  // fill rbuf with size and vector/array values

  rbuf[0] = n1;
  rbuf[1] = n2;
  if (vecflag) memcpy(&rbuf[2],vstore,n1*sizeof(double));
  else if (arrayflag) memcpy(&rbuf[2],&astore[0][0],n1*n2*sizeof(double));

  int n = n1*n2 + 2;
  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(rbuf,sizeof(double),n,fp);
  }
}

/* ----------------------------------------------------------------------
   use global array from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixStore::restart(char *buf)
{
  // first 2 values in buf are vec/array sizes

  double *dbuf = (double *) buf;
  int n1_restart = dbuf[0];
  int n2_restart = dbuf[1];

  // if size of vec/array has changed,
  //   means the restart file is setting size of vec or array and doing init
  //   because caller did not know size at time this fix was instantiated
  // reallocate vstore or astore accordingly

  if (n1 != n1_restart || n2 != n2_restart) {
    memory->destroy(vstore);
    memory->destroy(astore);
    memory->destroy(rbuf);
    vstore = nullptr;
    astore = nullptr;

    vecflag = arrayflag = 0;
    if (n2_restart == 1) vecflag = 1;
    else arrayflag = 1;
    n1 = n1_restart;
    n2 = n2_restart;
    if (vecflag) memory->create(vstore,n1,"fix/store:vstore");
    else if (arrayflag) memory->create(astore,n1,n2,"fix/store:astore");
    memory->create(rbuf,n1*n2+2,"fix/store:rbuf");
  }

  int n = n1*n2;
  if (vecflag) memcpy(vstore,&dbuf[2],n*sizeof(double));
  else if (arrayflag) memcpy(&astore[0][0],&dbuf[2],n*sizeof(double));
}

/* ----------------------------------------------------------------------
   allocate atom-based array
------------------------------------------------------------------------- */

void FixStore::grow_arrays(int nmax)
{
  if (vecflag) memory->grow(vstore,nmax,"store:vstore");
  else if (arrayflag) memory->grow(astore,nmax,n2,"store:astore");
  else if (tensorflag) memory->grow(tstore,nmax,n2,n3,"store:tstore");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */

void FixStore::copy_arrays(int i, int j, int /*delflag*/)
{
  if (disable) return;

  if (vecflag) {
    vstore[j] = vstore[i];
  } else if (arrayflag) {
    for (int m = 0; m < nvalues; m++)
      astore[j][m] = astore[i][m];
  } else if (tensorflag) {
    memcpy(&tstore[j][0][0],&tstore[i][0][0],nbytes);
  }
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixStore::pack_exchange(int i, double *buf)
{
  if (disable) return 0;

  if (vecflag) {
    buf[0] = vstore[i];
  } else if (arrayflag) {
    for (int m = 0; m < nvalues; m++)
      buf[m] = astore[i][m];
  } else if (tensorflag) {
    memcpy(buf,&tstore[i][0][0],nbytes);
  }

  return nvalues;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixStore::unpack_exchange(int nlocal, double *buf)
{
  if (disable) return 0;

  if (vecflag) {
    vstore[nlocal] = buf[0];
  } else if (arrayflag) {
    for (int m = 0; m < nvalues; m++)
      astore[nlocal][m] = buf[m];
  } else if (tensorflag) {
    memcpy(&tstore[nlocal][0][0],buf,nbytes);
  }

  return nvalues;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixStore::pack_restart(int i, double *buf)
{
  if (disable) {
    buf[0] = 0;
    return 1;
  }

  // pack buf[0] this way because other fixes unpack it
  buf[0] = nvalues+1;

  if (vecflag) {
    buf[1] = vstore[i];
  } else if (arrayflag) {
    for (int m = 0; m < nvalues; m++)
      buf[m+1] = astore[i][m];
  } else if (tensorflag) {
    memcpy(&buf[1],&tstore[i][0][0],nbytes);
  }

  return nvalues+1;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixStore::unpack_restart(int nlocal, int nth)
{
  if (disable) return;

  double **extra = atom->extra;

  // skip to Nth set of extra values
  // unpack the Nth first values this way because other fixes pack them

  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int> (extra[nlocal][m]);
  m++;

  if (vecflag) {
    vstore[nlocal] = extra[nlocal][m];
  } else if (arrayflag) {
    for (int i = 0; i < nvalues; i++)
      astore[nlocal][i] = extra[nlocal][m++];
  } else if (tensorflag) {
    memcpy(&tstore[nlocal][0][0],&extra[nlocal][m],nbytes);
  }
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixStore::maxsize_restart()
{
  if (disable) return 1;
  return nvalues+1;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixStore::size_restart(int /*nlocal*/)
{
  if (disable) return 1;
  return nvalues+1;
}

/* ----------------------------------------------------------------------
   memory usage of global or peratom atom-based array
------------------------------------------------------------------------- */

double FixStore::memory_usage()
{
  double bytes = 0.0;
  if (flavor == GLOBAL) bytes += n1*n2 * sizeof(double);
  if (flavor == PERATOM) bytes += atom->nmax*n2*n3 * sizeof(double);
  return bytes;
}
