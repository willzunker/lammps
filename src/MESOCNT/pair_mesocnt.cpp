#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include "pair_mesocnt.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "update.h"
#include "integrate.h"
#include "respa.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define MAXLINE 1024
#define SMALL 1.0e-6

/* ---------------------------------------------------------------------- */

PairMesoCNT::PairMesoCNT(LAMMPS *lmp) : Pair(lmp)
{
  respa_enable = 0;
  writedata = 1;

  
}

/* ---------------------------------------------------------------------- */

PairMesoCNT::~PairMesoCNT()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(gamma_data);
    memory->destroy(uinf_data);
    memory->destroy(starth_usemi);
    memory->destroy(startzeta_phi);
    memory->destroy(delh_usemi);
    memory->destroy(delzeta_phi);

    memory->destroy(usemi_data);
    memory->destroy(phi_data);
    memory->destroy(gamma_coeff);
    memory->destroy(uinf_coeff);
    
    memory->destroy(usemi_coeff);
    memory->destroy(phi_coeff);
 
    memory->destroy(redlist);
    memory->destroy(chain);
    memory->destroy(nchain);
    memory->destroy(end);

    memory->destroy(p1);
    memory->destroy(p2);
    memory->destroy(param);
    memory->destroy(flocal);
    memory->destroy(basis);
  }
}

/* ---------------------------------------------------------------------- */

void PairMesoCNT::compute(int eflag, int vflag)
{
  int inum,i,i1,i2,j,j1,j2,jj,jj1,jj2,j1num,j2num,listmax,numred,inflag,n;
  int ind,m,mm,loc1,loc2,cid,cnum,cn,tag1,tag2,idprev,idnext;
  double w,sumw,sumwreq,evdwl;
  double *r1,*r2,*q1,*q2,*qe;

  int *ilist,*j1list,*j2list,*numneigh,**firstneigh;
 
  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  int **bondlist = neighbor->bondlist;
  int *tag = atom->tag;
  int *mol = atom->molecule;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
 
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  for(n = 0; n < nbondlist; n++){
    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    
    r1 = x[i1];
    r2 = x[i2];
    
    // reduce neighbors to common list
    
    if(i1 > nlocal-1){
      j1num = 0;
    }
    else{
      j1list = firstneigh[i1];
      j1num = numneigh[i1];
    }
    if(i2 > nlocal-1){
      j2num = 0;
    }
    else{
      j2list = firstneigh[i2];
      j2num = numneigh[i2];
    }

    listmax = j1num + j2num;
    if(listmax < 2) continue;
    if(listmax > redlist_size){
      redlist_size = 2 * listmax;
      memory->grow(redlist,redlist_size,"pair:redlist");
    }

    numred = 0;
    for(jj1 = 0; jj1 < j1num; jj1++){
      ind = j1list[jj1];
      if ((mol[ind] == mol[i1]) && (abs(tag[ind] - tag[i1]) < 5)) continue;
      redlist[numred++] = ind;
    }
    for(jj2 = 0; jj2 < j2num; jj2++){
      for(jj1 = 0; jj1 < j1num; jj1++){
        if(j1list[jj1] == j2list[jj2]){
	  inflag = 1;
	  break;
	}
      }
      if(inflag){
        inflag = 0;
	continue;
      }
      ind = j2list[jj2];
      if ((mol[ind] == mol[i2]) && (abs(tag[ind] - tag[i2]) < 5)) continue; 
      redlist[numred++] = ind;
    }
    
    if (numred < 2) continue;
    
    // insertion sort according to atom-id
    
    sort(redlist,numred);
 
    // split into connected chains
    
    cid = 0;
    cnum = 0;
    
    if(numred > chain_size){
      chain_size = 2 * numred;
      memory->destroy(chain);
      memory->create(chain,chain_size,chain_size,"pair:chain");
      memory->grow(nchain,chain_size,"pair:nchain");
    }
    
    for(jj = 0; jj < numred-1; jj++){
      j1 = redlist[jj];
      j2 = redlist[jj+1];
      chain[cid][cnum++] = j1;
      if(abs(tag[j1] - tag[j2]) != 1 || mol[j1] != mol[j2]){
        nchain[cid++] = cnum;
	cnum = 0;
      }
    }
    chain[cid][cnum++] = redlist[numred-1];
    nchain[cid++] = cnum;

    // check for ends
 
    if(cid > end_size){
      end_size = 2 * cid;
      memory->grow(end,end_size,"pair:end_size");
    }
    
    for(i = 0; i < cid; i++){
      cn = nchain[i];
      loc1 = chain[i][0];
      loc2 = chain[i][cn-1];
      tag1 = tag[loc1];
      tag2 = tag[loc2];
      end[i] = 0;
      if(tag1 == 1) end[i] = 1;
      else{
        idprev = atom->map(tag1-1);
	if(idprev == -1 || mol[loc1] != mol[idprev]) end[i] = 1;
      }
      if(tag2 == atom->natoms) end[i] = 2;
      else{
        idnext = atom->map(tag2+1);
	if(idnext == -1 || mol[loc2] != mol[idnext]) end[i] = 2;
      }
    }

    // compute subsitute chains, forces and energies

    using namespace MathExtra;
    
    for(i = 0; i < cid; i++){
      if(nchain[i] < 2) continue;

      zero3(p1);
      zero3(p2);
      sumw = 0;

      loc1 = chain[i][0];
      loc2 = chain[i][nchain[i]-2];
      if(end[i] == 1) qe = x[loc1];
      if(end[i] == 2) qe = x[loc2];

      for(j = 0; j < nchain[i]-1; j++){
	q1 = x[chain[i][j]];
	q2 = x[chain[i][j+1]];

	w = weight(r1,r2,q1,q2);
	if(w < SMALL){ 
          if(end[i] == 1 && j == 0) end[i] = 0;
	  else if(end[i] == 2 && j == nchain[i] - 2) end[i] = 0;
	  continue;
	}
	sumw += w;
	MathExtra::scaleadd3(w,q1,p1,p1);
	MathExtra::scaleadd3(w,q2,p2,p2);
      }
      if (sumw < SMALL) continue;

      sumwreq = 1 / sumw;
      scale3(sumwreq,p1);
      scale3(sumwreq,p2);

      if(end[i] == 0){
        geominf(r1,r2,p1,p2,param,basis);
        if(param[0] > cutoff) continue;
        finf(param,flocal);
      }
      else{
        geomsemi(r1,r2,p1,p2,qe,param,basis);
	if(param[0] > cutoff) continue;
	fsemi(param,flocal);
      }
            
      f[i1][0] += flocal[0][0]*basis[0][0] 
	      + flocal[0][1]*basis[1][0] 
	      + flocal[0][2]*basis[2][0];
      f[i1][1] += flocal[0][0]*basis[0][1]
	      + flocal[0][1]*basis[1][1]
	      + flocal[0][2]*basis[2][1];
      f[i1][2] += flocal[0][0]*basis[0][2]
	      + flocal[0][1]*basis[1][2]
	      + flocal[0][2]*basis[2][2];
      f[i2][0] += flocal[1][0]*basis[0][0] 
	      + flocal[1][1]*basis[1][0] 
	      + flocal[1][2]*basis[2][0];
      f[i2][1] += flocal[1][0]*basis[0][1]
	      + flocal[1][1]*basis[1][1]
	      + flocal[1][2]*basis[2][1];
      f[i2][2] += flocal[1][0]*basis[0][2]
	      + flocal[1][1]*basis[1][2]
	      + flocal[1][2]*basis[2][2];
    
      if(eflag_either){
	if(eflag_global){
          if(end[i] == 0) evdwl = uinf(param);
          else evdwl = usemi(param);
	  eng_vdwl += 0.5 * evdwl;
	}
	if(eflag_atom){
	  eatom[i1] += 0.25 * evdwl;
	  eatom[i2] += 0.25 * evdwl;
	}
      }

      if(vflag_atom){
        vatom[i1][0] += f[i1][0]*x[i1][0];
	vatom[i1][1] += f[i1][1]*x[i1][1];
	vatom[i1][2] += f[i1][2]*x[i1][2];
	vatom[i1][3] += f[i1][1]*x[i1][0];
	vatom[i1][4] += f[i1][2]*x[i1][0];
	vatom[i1][5] += f[i1][2]*x[i1][1];
        vatom[i2][0] += f[i2][0]*x[i2][0];
	vatom[i2][1] += f[i2][1]*x[i2][1];
	vatom[i2][2] += f[i2][2]*x[i2][2];
	vatom[i2][3] += f[i2][1]*x[i2][0];
	vatom[i2][4] += f[i2][2]*x[i2][0];
	vatom[i2][5] += f[i2][2]*x[i2][1];
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

void PairMesoCNT::init_style()
{
  int irequest;
  irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

/* ---------------------------------------------------------------------- */

double PairMesoCNT::init_one(int i, int j)
{
  return cutoff;
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairMesoCNT::allocate()
{
  allocated = 1;
  int n = atom->ntypes;
  redlist_size = 64;
  chain_size = 64;
  end_size = 64;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");
 
  memory->create(gamma_data,gamma_points,"pair:gamma_data");
  memory->create(uinf_data,pot_points,"pair:uinf_data");
  memory->create(starth_usemi,pot_points,"pair:starth_usemi");
  memory->create(startzeta_phi,pot_points,"pair:startzeta_phi");
  memory->create(delh_usemi,pot_points,"pair:delh_usemi");
  memory->create(delzeta_phi,pot_points,"pair:delzeta_phi");
  
  memory->create(usemi_data,pot_points,pot_points,"pair:usemi_data");
  memory->create(phi_data,pot_points,pot_points,"pair:phi_data");
  memory->create(gamma_coeff,gamma_points,4,"pair:gamma_coeff");
  memory->create(uinf_coeff,pot_points,4,"pair:uinf_coeff");
  
  memory->create(usemi_coeff,pot_points,pot_points,4,"pair:usemi_coeff");
  memory->create(phi_coeff,pot_points,pot_points,4,"pair:phi_coeff");

  memory->create(redlist,redlist_size,"pair:redlist");
  memory->create(chain,chain_size,chain_size,"pair:chain");
  memory->create(nchain,chain_size,"pair:nchain");
  memory->create(end,end_size,"pair:end");

  memory->create(p1,3,"pair:p1");
  memory->create(p2,3,"pair:p2");
  memory->create(param,5,"pair:param");
  memory->create(flocal,2,3,"pair:flocal");
  memory->create(basis,3,3,"pair:basis");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairMesoCNT::settings(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR,"Illegal pair_style command");

  gamma_points = force->inumeric(FLERR,arg[0]);
  pot_points = force->inumeric(FLERR,arg[1]);
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairMesoCNT::coeff(int narg, char **arg)
{
  if (narg != 10) error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  n = force->inumeric(FLERR,arg[2]);
  sigma = force->numeric(FLERR,arg[3]);
  epsilon = force->numeric(FLERR,arg[4]);
  n_sigma = force->numeric(FLERR,arg[5]);

  gamma_file = arg[6];
  uinf_file = arg[7];
  usemi_file = arg[8];
  phi_file = arg[9];

  angstrom = force->angstrom;
  angstromrec = 1 / angstrom;
  qelectron = force->qelectron;
  qelectronrec = 1 / qelectron;
  forceunit = qelectron * angstromrec;

  radius = 1.421*3*n / MY_2PI * angstrom;
  radiussq = radius * radius;
  diameter = 2 * radius;
  rc = 3 * sigma;
  comega = 0.275*(1.0 - 1.0/(1.0 + 0.59*radius*angstromrec));
  ctheta = 0.35 + 0.0226*(radius*angstromrec - 6.785);

  //Parse and bcast data
  MPI_Comm_rank(world,&me);
  if (me == 0) { 
    read_file(gamma_file,gamma_data,&start_gamma,&del_gamma,gamma_points);
    read_file(uinf_file,uinf_data,&start_uinf,&del_uinf,pot_points);
    read_file(usemi_file,usemi_data,starth_usemi,&startxi_usemi,
		    delh_usemi,&delxi_usemi,pot_points);
    read_file(phi_file,phi_data,startzeta_phi,&starth_phi,
		    delzeta_phi,&delh_phi,pot_points);
  }
  
  MPI_Bcast(&start_gamma,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&del_gamma,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&start_uinf,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&del_uinf,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&startxi_usemi,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&delxi_usemi,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&starth_phi,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&delh_phi,1,MPI_DOUBLE,0,world);

  MPI_Bcast(gamma_data,gamma_points,MPI_DOUBLE,0,world);
  MPI_Bcast(uinf_data,pot_points,MPI_DOUBLE,0,world);
  MPI_Bcast(starth_usemi,pot_points,MPI_DOUBLE,0,world);
  MPI_Bcast(startzeta_phi,pot_points,MPI_DOUBLE,0,world); 
  MPI_Bcast(delh_usemi,pot_points,MPI_DOUBLE,0,world);
  MPI_Bcast(delzeta_phi,pot_points,MPI_DOUBLE,0,world);
  for(int i = 0; i < pot_points; i++){
    MPI_Bcast(usemi_data[i],pot_points,MPI_DOUBLE,0,world);
    MPI_Bcast(phi_data[i],pot_points,MPI_DOUBLE,0,world);
  }

  spline_coeff(gamma_data,gamma_coeff,gamma_points);
  spline_coeff(uinf_data,uinf_coeff,pot_points);
  spline_coeff(usemi_data,usemi_coeff,pot_points);
  spline_coeff(phi_data,phi_coeff,pot_points);

  int n = atom->ntypes;

  cutoff = rc + diameter;
  cutoffsq = cutoff * cutoff;
  
  for (int i = 1; i <= n; i++){
    for (int j = i; j <= n; j++){
      setflag[i][j] = 1;
      cutsq[i][j] = cutoffsq;
    }
  }
}

/* ---------------------------------------------------------------------- */

double PairMesoCNT::spline(double x, double xstart, double dx, 
		double **coeff, int coeff_size)
{
  int i = floor((x - xstart)/dx); 
  if(i < 0){
    i = 0;
    x = xstart;
  }
  else if(i > coeff_size-2){ 
    i = coeff_size - 2;
    x = xstart + (coeff_size - 1)*dx;
  }
 
  double xlo = xstart + i*dx;
  double xbar = (x - xlo)/dx;

  return coeff[i][0] + xbar*(coeff[i][1] 
		  + xbar*(coeff[i][2] + xbar*coeff[i][3]));
}

/* ---------------------------------------------------------------------- */

double PairMesoCNT::spline(double x, double y, double *xstart, double ystart,
		double *dx, double dy, double ***coeff, int coeff_size)
{
  int i = floor((y - ystart)/dy);
  if(i < 0){
    i = 0;
    y = ystart;
  }
  else if(i > coeff_size-3){
    i = coeff_size - 3;
    y = ystart + (coeff_size - 2)*dy;
  }
  
  double ylo = ystart + i*dy;
  double ybar = (y - ylo)/dy;

  // compute coefficients in y
  
  double a0, a1, a2, a3;
  double p0, p1, p2, p3;
  
  p1 = spline(x,xstart[i],dx[i],coeff[i],coeff_size);
  p2 = spline(x,xstart[i+1],dx[i+1],coeff[i+1],coeff_size);
  
  a0 = p1;

  if(i == 0){
    p3 = spline(x,xstart[i+2],dx[i+2],coeff[i+2],coeff_size);

    a1 = p2 - p1;
    a3 = 0.5*(p3 - 2*p2 + p1);
    a2 = -a3;
  }
  else if(i == coeff_size-2){
    p0 = spline(x,xstart[i-1],dx[i-1],coeff[i-1],coeff_size);
    
    a1 = 0.5*(p2 - p0);
    a3 = 0.5*(p2 - 2*p1 + p0);
    a2 = -2*a3;
  }
  else{
    p0 = spline(x,xstart[i-1],dx[i-1],coeff[i-1],coeff_size);
    p3 = spline(x,xstart[i+2],dx[i+2],coeff[i+2],coeff_size);

    a1 = 0.5*(p2 - p0);
    a2 = 0.5*(-p3 + 4*p2 - 5*p1 + 2*p0);
    a3 = 0.5*(p3 - 3*p2 + 3*p1 - p0);
  }

  return a0 + ybar*(a1 + ybar*(a2 + a3*ybar));
}

/* ---------------------------------------------------------------------- */

double PairMesoCNT::dspline(double x, double xstart, double dx, 
		double **coeff, int coeff_size)
{
  int i = floor((x - xstart)/dx);
  if(i < 0){
    i = 0;
    x = xstart;
  }
  else if(i > coeff_size-2){
    i = coeff_size - 2;
    x = xstart + (coeff_size - 1)*dx;
  }
 
  double xlo = xstart + i*dx;
  double xbar = (x - xlo)/dx;

  return (coeff[i][1] + xbar*(2*coeff[i][2] + 3*xbar*coeff[i][3])) / dx;
}

/* ---------------------------------------------------------------------- */

double PairMesoCNT::dxspline(double x, double y, double *xstart, double ystart,
		double *dx, double dy, double ***coeff, int coeff_size)
{
  int i = floor((y - ystart)/dy);
  if(i < 0){
    i = 0;
    y = ystart;
  }
  if(i > coeff_size-3){ 
    i = coeff_size - 3;
    y = ystart + (coeff_size - 2)*dy;
  }

  double ylo = ystart + i*dy;
  double ybar = (y - ylo)/dy;

  // compute coefficients in y
  
  double a0, a1, a2, a3;
  double p0, p1, p2, p3;
  
  p1 = dspline(x,xstart[i],dx[i],coeff[i],coeff_size);
  p2 = dspline(x,xstart[i+1],dx[i+1],coeff[i+1],coeff_size);
  
  a0 = p1;

  if(i == 0){
    p3 = dspline(x,xstart[i+2],dx[i+2],coeff[i+2],coeff_size);

    a1 = p2 - p1;
    a3 = 0.5*(p3 - 2*p2 + p1);
    a2 = -a3;
  }
  else if(i == coeff_size-2){
    p0 = dspline(x,xstart[i-1],dx[i-1],coeff[i-1],coeff_size);
    
    a1 = 0.5*(p2 - p0);
    a3 = 0.5*(p2 - 2*p1 + p0);
    a2 = -2*a3;
  }
  else{
    p0 = dspline(x,xstart[i-1],dx[i-1],coeff[i-1],coeff_size);
    p3 = dspline(x,xstart[i+2],dx[i+2],coeff[i+2],coeff_size);

    a1 = 0.5*(p2 - p0);
    a2 = 0.5*(-p3 + 4*p2 - 5*p1 + 2*p0);
    a3 = 0.5*(p3 - 3*p2 + 3*p1 - p0);
  }

  return a0 + ybar*(a1 + ybar*(a2 + a3*ybar));
}

/* ---------------------------------------------------------------------- */

double PairMesoCNT::dyspline(double x, double y, double *xstart, double ystart,
		double *dx, double dy, double ***coeff, int coeff_size)
{
  int i = floor((y - ystart)/dy);
  if(i < 0){
    i = 0;
    y = ystart;
  }
  if(i > coeff_size-3){
    i = coeff_size - 3;
    y = ystart + (coeff_size - 2)*dy;
  }
 
  double ylo = ystart + i*dy;
  double ybar = (y - ylo)/dy;

  // compute coefficients in y
  
  double a0, a1, a2, a3;
  double p0, p1, p2, p3;
  
  p1 = spline(x,xstart[i],dx[i],coeff[i],coeff_size);
  p2 = spline(x,xstart[i+1],dx[i+1],coeff[i+1],coeff_size);
  
  a0 = p1;

  if(i == 0){
    p3 = spline(x,xstart[i+2],dx[i+2],coeff[i+2],coeff_size);

    a1 = p2 - p1;
    a3 = 0.5*(p3 - 2*p2 + p1);
    a2 = -a3;
  }
  else if(i == coeff_size-2){
    p0 = spline(x,xstart[i-1],dx[i-1],coeff[i-1],coeff_size);
    
    a1 = 0.5*(p2 - p0);
    a3 = 0.5*(p2 - 2*p1 + p0);
    a2 = -2*a3;
  }
  else{
    p0 = spline(x,xstart[i-1],dx[i-1],coeff[i-1],coeff_size);
    p3 = spline(x,xstart[i+2],dx[i+2],coeff[i+2],coeff_size);

    a1 = 0.5*(p2 - p0);
    a2 = 0.5*(-p3 + 4*p2 - 5*p1 + 2*p0);
    a3 = 0.5*(p3 - 3*p2 + 3*p1 - p0);
  }

  return (a1 + ybar*(2*a2 + 3*a3*ybar)) / dy;
}

/* ---------------------------------------------------------------------- */

void PairMesoCNT::spline_coeff(double *data, double **coeff, int data_size)
{
  for(int i = 0; i < data_size-1; i++){
    if(i == 0){
      coeff[i][0] = data[i];
      coeff[i][1] = data[i+1] - data[i];
      coeff[i][3] = 0.5*(data[i+2] - 2*data[i+1] + data[i]);
      coeff[i][2] = -coeff[i][3];
    }
    else if(i == data_size-2){
      coeff[i][0] = data[i];
      coeff[i][1] = 0.5*(data[i+1] - data[i-1]);
      coeff[i][3] = 0.5*(-data[i+1] + 2*data[i] - data[i-1]);
      coeff[i][2] = -2*coeff[i][3];
    }
    else{
      coeff[i][0] = data[i];
      coeff[i][1] = 0.5*(data[i+1] - data[i-1]);
      coeff[i][2] = 0.5*(-data[i+2] + 4*data[i+1] - 5*data[i] + 2*data[i-1]);
      coeff[i][3] = 0.5*(data[i+2] - 3*data[i+1] + 3*data[i] - data[i-1]);
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairMesoCNT::spline_coeff(double **data, double ***coeff, int data_size)
{
  for(int i = 0; i < data_size; i++){
    for(int j = 0; j < data_size-1; j++){
      if(j == 0){
        coeff[i][j][0] = data[j][i];
        coeff[i][j][1] = data[j+1][i] - data[j][i];
        coeff[i][j][3] = 0.5*(data[j+2][i] - 2*data[j+1][i] + data[j][i]);
        coeff[i][j][2] = -coeff[i][j][3];
      }
      else if(j == data_size-2){
        coeff[i][j][0] = data[j][i];
        coeff[i][j][1] = 0.5*(data[j+1][i] - data[j-1][i]);
        coeff[i][j][3] = 0.5*(-data[j+1][i] + 2*data[j][i] - data[j-1][i]);
        coeff[i][j][2] = -2*coeff[i][j][3];
      }
      else{
        coeff[i][j][0] = data[j][i];
        coeff[i][j][1] = 0.5*(data[j+1][i] - data[j-1][i]);
        coeff[i][j][2] = 0.5*(-data[j+2][i] + 4*data[j+1][i] 
			- 5*data[j][i] + 2*data[j-1][i]);
        coeff[i][j][3] = 0.5*(data[j+2][i] - 3*data[j+1][i] 
			+ 3*data[j][i] - data[j-1][i]);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairMesoCNT::read_file(char *file, double *data, 
		double *startx, double *dx, int ninput)
{
  char line[MAXLINE];

  // open file

  FILE *fp = force->open_potential(file);
  if (fp == NULL) {
    std::string str("Cannot open file ");
    str += file;
    error->one(FLERR,str.c_str());
  }

  // read values from file
  
  int cerror = 0;
  int serror = 0;
  double x,xtemp;

  for(int i = 0; i < ninput; i++){
    if(i > 0) xtemp = x;
    if(NULL == fgets(line,MAXLINE,fp)){
      std::string str("Premature end of file in pair table ");
      str += file;
      error->one(FLERR,str.c_str());
    }
    if(2 != sscanf(line,"%lg %lg",&x, &data[i])) ++cerror; 
    if(i == 0) *startx = x;
    if(i > 0){
      if(i == 1) *dx = x - xtemp;
      if((*dx -  x + xtemp) / *dx > SMALL) ++serror;
    }
  }

  // warn if data was read incompletely, e.g. columns were missing

  if (cerror) {
    char str[128];
    sprintf(str,"%d of %d lines in table were incomplete\n"
            "  or could not be parsed completely",cerror,ninput);
    error->warning(FLERR,str);
  }

  // warn if spacing between data points is not constant
  
  if (serror) {
    char str[128];
    sprintf(str, "%d spacings were different\n"
	    "  from first entry",serror);
    error->warning(FLERR,str);
  }

}

/* ---------------------------------------------------------------------- */

void PairMesoCNT::read_file(char *file, double **data, 
		double *startx, double *starty, double *dx, double *dy,
		int ninput)
{
  char line[MAXLINE];

  // open file

  FILE *fp = force->open_potential(file);
  if (fp == NULL) {
    std::string str("Cannot open file ");
    str += file;
    error->one(FLERR,str.c_str());
  }

  // read values from file
  
  int cerror = 0;
  int sxerror = 0;
  int syerror = 0;
  double x,y,xtemp,ytemp;

  for(int i = 0; i < ninput; i++){ 
    if(i > 0) ytemp = y;
    for(int j = 0; j < ninput; j++){
      if(j > 0) xtemp = x;
      if(NULL == fgets(line,MAXLINE,fp)){
	std::string str("Premature end of file in pair table ");
	str += file;
        error->one(FLERR,str.c_str());
      }
      if(3 != sscanf(line,"%lg %lg %lg",&x,&y,&data[j][i])) ++cerror; 
      if(j == 0) startx[i] = x;
      if(j > 0){
        if(j == 1) dx[i] = x - xtemp;
        if((dx[i] - x + xtemp)/dx[i] > SMALL) ++sxerror;
      }
    }
    if(i == 0) *starty = y;
    if(i > 0){
      if(i == 1) *dy = y - ytemp;
      if((*dy - y + ytemp)/ *dy > SMALL) ++syerror;
    }
  }

  // warn if data was read incompletely, e.g. columns were missing

  if (cerror) {
    char str[128];
    sprintf(str,"%d of %d lines in table were incomplete\n"
            "  or could not be parsed completely",cerror,ninput);
    error->warning(FLERR,str);
  }
  
  // warn if spacing between data points is not constant
  
  if (sxerror) {
    char str[128];
    sprintf(str, "%d spacings in first column were different\n"
	    "  from first block entries",sxerror);
    error->warning(FLERR,str);
  }

  if (syerror) {
    char str[128];
    sprintf(str, "%d spacings in second column were different\n"
	    "  from first entry",syerror);
    error->warning(FLERR,str);
  }

}

/* ---------------------------------------------------------------------- */

double PairMesoCNT::uinf(double *param)
{
  double h = param[0] * angstromrec;
  double alpha = param[1];
  double xi1 = param[2] * angstromrec;
  double xi2 = param[3] * angstromrec;

  double sin_alpha = sin(alpha);
  double sin_alphasq = sin_alpha*sin_alpha;
  if(sin_alphasq < SMALL){
    return (xi2 - xi1) * spline(h,start_uinf,del_uinf,uinf_coeff,pot_points)
	    * qelectron;
  }
  else{
    double omega = 1.0 / (1.0 - comega*sin_alphasq);
    double a = omega * sin_alpha;
    double zeta1 = xi1 * a;
    double zeta2 = xi2 * a;

    double phi1, phi2;
    if(zeta1 < 0) phi1 = -spline(-zeta1,h,startzeta_phi,starth_phi,
		    delzeta_phi,delh_phi,phi_coeff,pot_points);
    else phi1 = spline(zeta1,h,startzeta_phi,starth_phi,
		    delzeta_phi,delh_phi,phi_coeff,pot_points);

    if(zeta2 < 0) phi2 = -spline(-zeta2,h,startzeta_phi,starth_phi,
		    delzeta_phi,delh_phi,phi_coeff,pot_points);
    else phi2 = spline(zeta2,h,startzeta_phi,starth_phi,
		    delzeta_phi,delh_phi,phi_coeff,pot_points);

    double gamma_orth = spline(h,start_gamma,del_gamma,
		    gamma_coeff,gamma_points);
    double gamma = 1.0 + (gamma_orth - 1.0)*sin_alphasq;

    return gamma * (phi2 - phi1) * qelectron / a;
  }
}

/* ---------------------------------------------------------------------- */

double PairMesoCNT::usemi(double *param)
{
  double h = param[0] * angstromrec;
  double alpha = param[1];
  double xi1 = param[2] * angstromrec;
  double xi2 = param[3] * angstromrec;
  double etae = param[4] * angstromrec;

  double sin_alpha = sin(alpha);
  double sin_alphasq = sin_alpha*sin_alpha;
  double cos_alpha = cos(alpha);
  double omega = 1.0 / (1.0 - comega*sin_alphasq);
  double theta = 1.0 - ctheta*sin_alphasq; 

  double a1 = omega * sin_alpha;
  double a2 = theta * etae;

  int points = 100;
  double delxi = (xi2 - xi1) / (points - 1);
  
  double sum = 0;

  for(int i = 0; i < points; i++){
    double xibar = xi1 + i*delxi;
    double g = xibar * a1;
    double hbar = sqrt(h*h + g*g);
    double zetabar = xibar*cos_alpha - a2;
    
    double c = 1.0;
    if(i == 0 || i == points-1) c = 0.5;

    sum += c * spline(hbar,zetabar,starth_usemi,startxi_usemi,
		  delh_usemi,delxi_usemi,usemi_coeff,pot_points);
  }

  double gamma_orth = spline(h,start_gamma,del_gamma,
		    gamma_coeff,gamma_points);
  double gamma = 1.0 + (gamma_orth - 1.0)*sin_alphasq;

  return delxi * gamma * sum * qelectron;
}

/* ---------------------------------------------------------------------- */

void PairMesoCNT::finf(double *param, double **f)
{
  double h = param[0] * angstromrec;
  double alpha = param[1];
  double xi1 = param[2] * angstromrec;
  double xi2 = param[3] * angstromrec;

  double sin_alpha = sin(alpha);
  double sin_alphasq = sin_alpha*sin_alpha;
  
  if(sin_alphasq < SMALL){
    f[0][0] = 0.5 * (xi2 - xi1) * dspline(h,start_uinf,del_uinf,
        uinf_coeff,pot_points) * forceunit;
    f[1][0] = f[0][0];
    f[0][1] = 0;
    f[1][1] = 0;
    f[0][2] = spline(h,start_uinf,del_uinf,uinf_coeff,pot_points)
	    * forceunit;
    f[1][2] = -f[0][2];
    return;
  }
  else{
    double sin_alpharec = 1.0 / sin_alpha;
    double sin_alpharecsq = sin_alpharec * sin_alpharec;
    double cos_alpha = cos(alpha);
    double cot_alpha = cos_alpha * sin_alpharec;

    double omega = 1.0 / (1.0 - comega*sin_alphasq);
    double a1 = omega * sin_alpha;
    double a1rec = 1.0 / a1;
    double domega = 2 * comega * cos_alpha * a1 * omega;
    
    double gamma_orth = spline(h,start_gamma,del_gamma,
		    gamma_coeff,gamma_points);
    double gamma = 1.0 + (gamma_orth - 1.0)*sin_alphasq;
    double gammarec = 1.0 / gamma;
    double dalpha_gamma = 2 * (gamma_orth - 1) * sin_alpha * cos_alpha;
    double dh_gamma = dspline(h,start_gamma,del_gamma,
		    gamma_coeff,gamma_points) * sin_alphasq;

    double zeta1 = xi1 * a1;
    double zeta2 = xi2 * a1;
    double hsq = h * h;
    double cutoffsq_local = cutoffsq * angstromrec * angstromrec;
    double phi1,phi2,dzeta_phi1,dzeta_phi2,dh_phi1,dh_phi2;
    
    if(zeta1 < 0){
      phi1 = -spline(-zeta1,h,startzeta_phi,starth_phi,
      	delzeta_phi,delh_phi,phi_coeff,pot_points);
      dzeta_phi1 = dxspline(-zeta1,h,startzeta_phi,starth_phi,
	delzeta_phi,delh_phi,phi_coeff,pot_points);
      dh_phi1 = -dyspline(-zeta1,h,startzeta_phi,starth_phi,
  	delzeta_phi,delh_phi,phi_coeff,pot_points);
    }
    else{ 
      phi1 = spline(zeta1,h,startzeta_phi,starth_phi,
	delzeta_phi,delh_phi,phi_coeff,pot_points);
      dzeta_phi1 = dxspline(zeta1,h,startzeta_phi,starth_phi,
	delzeta_phi,delh_phi,phi_coeff,pot_points);
      dh_phi1 = dyspline(zeta1,h,startzeta_phi,starth_phi,
  	delzeta_phi,delh_phi,phi_coeff,pot_points);
    }
    if(zeta2 < 0){
      phi2 = -spline(-zeta2,h,startzeta_phi,starth_phi,
        delzeta_phi,delh_phi,phi_coeff,pot_points);
      dzeta_phi2 = dxspline(-zeta2,h,startzeta_phi,starth_phi,
	delzeta_phi,delh_phi,phi_coeff,pot_points);
      dh_phi2 = -dyspline(-zeta2,h,startzeta_phi,starth_phi,
  	delzeta_phi,delh_phi,phi_coeff,pot_points);
    }
    else{ 
      phi2 = spline(zeta2,h,startzeta_phi,starth_phi,
	delzeta_phi,delh_phi,phi_coeff,pot_points);
      dzeta_phi2 = dxspline(zeta2,h,startzeta_phi,starth_phi,
	delzeta_phi,delh_phi,phi_coeff,pot_points);
      dh_phi2 = dyspline(zeta2,h,startzeta_phi,starth_phi,
  	delzeta_phi,delh_phi,phi_coeff,pot_points);
    }

    double diff_dzeta_phi = dzeta_phi2 - dzeta_phi1;
    
    double a2 = gamma * a1rec;
    double u = a2 * (phi2 - phi1);
    double a3 = u * gammarec;
  
    double dh_u = dh_gamma * a3 + a2 * (dh_phi2 - dh_phi1);
    double dalpha_u = dalpha_gamma * a3
	  + a1rec * (domega*sin_alpha + omega*cos_alpha)
	  * (gamma*(xi2*dzeta_phi2 - xi1*dzeta_phi1) - u);

    double lrec = 1.0 / (xi2 - xi1);
    double cx = h * gamma * sin_alpharecsq * diff_dzeta_phi;
    double cy = gamma * cot_alpha * diff_dzeta_phi;

    f[0][0] = lrec * (xi2*dh_u - cx) * forceunit;
    f[1][0] = lrec * (-xi1*dh_u + cx) * forceunit;
    f[0][1] = lrec * (dalpha_u - xi2*cy) * forceunit;
    f[1][1] = lrec * (-dalpha_u + xi1*cy) * forceunit;
    f[0][2] = gamma * dzeta_phi1 * forceunit;
    f[1][2] = -gamma * dzeta_phi2 * forceunit;
  }
}

/* ---------------------------------------------------------------------- */

void PairMesoCNT::fsemi(double *param, double **f)
{
  double h = param[0] * angstromrec;
  double alpha = param[1];
  double xi1 = param[2] * angstromrec;
  double xi2 = param[3] * angstromrec;
  double etae = param[4] * angstromrec;

  double sin_alpha = sin(alpha);
  double sin_alphasq = sin_alpha*sin_alpha;
  double cos_alpha = cos(alpha);

  double omega = 1.0 / (1.0 - comega*sin_alphasq);
  double omegasq = omega * omega;
  double domega = 2 * comega * sin_alpha * cos_alpha * omegasq;

  double theta = 1.0 - ctheta*sin_alphasq; 
  double dtheta = -2 * ctheta * sin_alpha * cos_alpha;

  double a1 = omega * sin_alpha;
  double a1sq = a1 * a1;
  double a2 = theta * etae;

  double gamma_orth = spline(h,start_gamma,del_gamma,
		    gamma_coeff,gamma_points);
  double gamma = 1.0 + (gamma_orth - 1.0)*sin_alphasq;
  double gammarec = 1.0 / gamma;
  double dalpha_gamma = 2 * (gamma_orth - 1) * sin_alpha * cos_alpha;
  double dh_gamma = dspline(h,start_gamma,del_gamma,
		    gamma_coeff,gamma_points) * sin_alphasq;
 
  int points = 100;
  double delxi = (xi2 - xi1) / (points - 1);
  double a3 = delxi * gamma;
  
  double jh = 0;
  double jh1 = 0;
  double jh2 = 0;
  double jxi = 0;
  double jxi1 = 0;
  double ubar = 0;

  for(int i = 0; i < points; i++){
    double xibar = xi1 + i*delxi;
    double g = xibar * a1;
    double hbar = sqrt(h*h + g*g);
    double zetabar = xibar*cos_alpha - a2;

    double c = 1.0;
    if(i == 0 || i == points-1) c = 0.5;

    double u = c * spline(hbar,zetabar,starth_usemi,startxi_usemi,
	delh_usemi,delxi_usemi,usemi_coeff,pot_points);
    double uh = c / hbar * dxspline(hbar,zetabar,starth_usemi,startxi_usemi,
	delh_usemi,delxi_usemi,usemi_coeff,pot_points);
    double uxi = c * dyspline(hbar,zetabar,starth_usemi,startxi_usemi,
	delh_usemi,delxi_usemi,usemi_coeff,pot_points);

    double uh1 = xibar * uh;
    jh += uh;
    jh1 += uh1;
    jh2 += xibar * uh1;
    jxi += uxi;
    jxi1 += xibar * uxi;
    ubar += u;
  }

  jh *= a3;
  jh1 *= a3;
  jh2 *= a3;
  jxi *= a3;
  jxi1 *= a3;
  ubar *= a3;

  double a4 = gammarec * ubar;
  double dh_ubar = dh_gamma*a4 + h*jh;
  double dalpha_ubar = dalpha_gamma*a4 
    + a1*(domega*sin_alpha + omega*cos_alpha)*jh2
    - sin_alpha * jxi1 - dtheta*etae*jxi;

  double cx = h * (omegasq*jh1 + cos_alpha*ctheta*jxi);
  double cy = sin_alpha * (cos_alpha*omegasq*jh1 + (ctheta-1)*jxi);
  double cz1 = a1sq*jh1 + cos_alpha*jxi;
  double cz2 = a1sq*jh2 + cos_alpha*jxi1;

  double lrec = 1.0 / (xi2 - xi1);
  f[0][0] = lrec * (xi2*dh_ubar - cx) * forceunit;
  f[1][0] = lrec * (cx - xi1*dh_ubar) * forceunit;
  f[0][1] = lrec * (dalpha_ubar - xi2*cy) * forceunit;
  f[1][1] = lrec * (xi1*cy - dalpha_ubar) * forceunit;
  f[0][2] = lrec * (cz2 + ubar - xi2*cz1) * forceunit;
  f[1][2] = lrec * (xi1*cz1 - cz2 - ubar) * forceunit;
}

/* ---------------------------------------------------------------------- */

void PairMesoCNT::geominf(const double *r1, const double *r2, 
		const double *p1, const double *p2, 
		double *param, double **basis)
{
  using namespace MathExtra;
  double r[3], p[3], delr[3], l[3], m[3], rbar[3], pbar[3], delrbar[3];
  double psil[3], psim[3], dell_psim[3], delpsil_m[3];
  double delr1[3], delr2[3], delp1[3], delp2[3];
  double *ex, *ey, *ez;
  double psi, denom, frac, taur, taup;
  double h, alpha, xi1, xi2;

  ex = basis[0];
  ey = basis[1];
  ez = basis[2];

  add3(r1,r2,r);
  scale3(0.5,r);
  add3(p1,p2,p);
  scale3(0.5,p);
  
  sub3(p,r,delr);

  sub3(r2,r1,l);
  normalize3(l,l);
  sub3(p2,p1,m);
  normalize3(m,m);

  psi = dot3(l,m);
  if(psi > 1.0) psi = 1.0;
  else if(psi < -1.0) psi = -1.0;
  denom = 1.0 - psi*psi;

  copy3(l,psil);
  scale3(psi,psil);
  copy3(m,psim);
  scale3(psi,psim);

  if(denom < SMALL){
    taur = dot3(delr,l);
    taup = 0;
  }
  else{
    frac = 1.0 / denom;
    sub3(l,psim,dell_psim);
    sub3(psil,m,delpsil_m);
    taur = dot3(delr,dell_psim) * frac;
    taup = dot3(delr,delpsil_m) * frac;
  }

  scaleadd3(taur,l,r,rbar);
  scaleadd3(taup,m,p,pbar);
  sub3(pbar,rbar,delrbar);
  
  h = len3(delrbar);
    
  copy3(delrbar,ex);
  copy3(l,ez);
  scale3(1/h,ex);
  cross3(ez,ex,ey);

  if(dot3(m,ey) < 0) alpha = acos(psi);
  else alpha = MY_2PI - acos(psi);

  sub3(r1,rbar,delr1);
  sub3(r2,rbar,delr2);
  xi1 = dot3(delr1,l);
  xi2 = dot3(delr2,l);

  param[0] = h;
  param[1] = alpha;
  param[2] = xi1;
  param[3] = xi2;
}

/* ---------------------------------------------------------------------- */

void PairMesoCNT::geomsemi(const double *r1, const double *r2, 
		const double *p1, const double *p2, const double *qe,
		double *param, double **basis)
{
  using namespace MathExtra;
  double r[3], p[3], delr[3], l[3], m[3], rbar[3], pbar[3], delrbar[3];
  double psil[3], psim[3], dell_psim[3], delpsil_m[3];
  double delr1[3], delr2[3], delp1[3], delp2[3], delpqe[3];
  double *ex, *ey, *ez;
  double psi, denom, frac, taur, taup, rhoe;
  double h, alpha, xi1, xi2, etae;

  ex = basis[0];
  ey = basis[1];
  ez = basis[2];

  add3(r1,r2,r);
  scale3(0.5,r);
  add3(p1,p2,p);
  scale3(0.5,p);
  
  sub3(p,r,delr);

  sub3(r2,r1,l);
  normalize3(l,l);
  sub3(p2,p1,m);
  normalize3(m,m);

  psi = dot3(l,m);
  if(psi > 1.0) psi = 1.0;
  else if(psi < -1.0) psi = -1.0;
  denom = 1.0 - psi*psi;

  copy3(l,psil);
  scale3(psi,psil);
  copy3(m,psim);
  scale3(psi,psim);

  sub3(p,qe,delpqe);
  rhoe = dot3(delpqe,m); 

  if(denom < SMALL){
    taur = dot3(delr,l) - rhoe*dot3(m,l);
    taup = -rhoe;
    etae = 0;
  }
  else{
    frac = 1.0 / denom;
    sub3(l,psim,dell_psim);
    sub3(psil,m,delpsil_m);
    taur = dot3(delr,dell_psim) * frac;
    taup = dot3(delr,delpsil_m) * frac;
    etae = -rhoe - taup;
  }

  scaleadd3(taur,l,r,rbar);
  scaleadd3(taup,m,p,pbar);
  sub3(pbar,rbar,delrbar);
  
  h = len3(delrbar);
    
  copy3(delrbar,ex);
  copy3(l,ez);
  scale3(1/h,ex);
  cross3(ez,ex,ey);

  if(dot3(m,ey) < 0) alpha = acos(psi);
  else alpha = MY_2PI - acos(psi);

  sub3(r1,rbar,delr1);
  sub3(r2,rbar,delr2);
  xi1 = dot3(delr1,l);
  xi2 = dot3(delr2,l);

  param[0] = h;
  param[1] = alpha;
  param[2] = xi1;
  param[3] = xi2;
  param[4] = etae;
}


/* ---------------------------------------------------------------------- */

double PairMesoCNT::weight(const double *r1, const double *r2, 
		const double *p1, const double *p2)
{
  using namespace MathExtra;
  double r[3], p[3];
  double rho, rhoc, rhomin;
  
  add3(r1,r2,r);
  add3(p1,p2,p);
  
  rhoc = sqrt(0.25*distsq3(r1,r2) + radiussq) 
	  + sqrt(0.25*distsq3(p1,p2) + radiussq) + rc;
  rhomin = 0.72 * rhoc;
  rho = 0.5 * sqrt(distsq3(r,p));

  return s((rho - rhomin)/(rhoc - rhomin));
}

/* ---------------------------------------------------------------------- */

int PairMesoCNT::heaviside(double x)
{
  if(x < 0) return 0;
  else return 1;
}

/* ---------------------------------------------------------------------- */

double PairMesoCNT::s(double x)
{
  return heaviside(-x) + heaviside(x)*heaviside(1-x)*(1 - x*x*(3 - 2*x));
}

/* ---------------------------------------------------------------------- */

void PairMesoCNT::sort(int *list, int size){
  int i,j,loc1,loc2;
  int *tag = atom->tag;
  for(i = 1; i < size; i++){
    j = i;
    loc1 = list[j-1];
    loc2 = list[j];
    while(j > 0 && tag[loc1] > tag[loc2]){
      list[j] = loc1;
      list[j-1] = loc2;
      j--;
      loc1 = list[j-1];
      loc2 = list[j];
    }
  } 
}
