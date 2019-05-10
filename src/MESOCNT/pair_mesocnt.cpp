#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
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
#include "utils.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define MAXLINE 1024

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
    memory->destroy(gamma_data);
    memory->destroy(uinf_data);
    memory->destroy(delh_usemi);
    memory->destroy(delzeta_phi);

    memory->destroy(usemi_data);
    memory->destroy(phi_data);
    memory->destroy(gamma_coeff);
    memory->destroy(uinf_coeff);
    
    memory->destroy(usemi_coeff);
    memory->destroy(phi_coeff);
  }
}

/* ---------------------------------------------------------------------- */

void PairMesoCNT::compute(int eflag, int vflag)
{

}

/* ---------------------------------------------------------------------- */



/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairMesoCNT::allocate()
{
  allocated = 1;
  
  memory->create(gamma_data,gamma_points,"pair:gamma_data");
  memory->create(uinf_data,pot_points,"pair:uinf_data");
  memory->create(delh_usemi,pot_points,"pair:delh_usemi");
  memory->create(delzeta_phi,pot_points,"pair:delzeta_phi");
  
  memory->create(usemi_data,pot_points,pot_points,"pair:usemi_data");
  memory->create(phi_data,pot_points,pot_points,"pair:phi_data");
  memory->create(gamma_coeff,gamma_points-1,4,"pair:gamma_coeff");
  memory->create(uinf_coeff,pot_points-1,4,"pair:uinf_coeff");
  
  memory->create(usemi_coeff,pot_points,pot_points-1,4,"pair:usemi_coeff");
  memory->create(phi_coeff,pot_points,pot_points-1,4,"pair:phi_coeff");
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

  radius = 1.421*3*n / MY_2PI;
  comega = 0.275*(1.0 - 1.0/(1.0 + 0.59*radius));
  ctheta = 0.35 + 0.0226*(radius - 6.785);

  //Parse and bcast data
  int me;
  MPI_Comm_rank(world,&me);
  if (me == 0) { 
    read_file(gamma_file,gamma_data,&del_gamma,gamma_points);
    read_file(uinf_file,uinf_data,&del_uinf,pot_points);
    read_file(usemi_file,usemi_data,delh_usemi,&delxi_usemi,pot_points);
    read_file(phi_file,phi_data,delzeta_phi,&delh_phi,pot_points);
  }
  
  MPI_Bcast(gamma_data,gamma_points,MPI_DOUBLE,0,world);
  MPI_Bcast(uinf_data,pot_points,MPI_DOUBLE,0,world); 
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
}

/* ---------------------------------------------------------------------- */

double PairMesoCNT::spline(double x, double xstart, double dx, 
		double **coeff, int coeff_size)
{
  int i = floor((x - xstart)/dx); 
  if(i < 0){
    i = 0;
    // warn if argument below spline range
    char str[128];
    sprintf(str,"Argument below spline interval");
    error->warning(FLERR,str);
  }
  else if(i > coeff_size-1){ 
    i = coeff_size-1;
    // warn if argument above spline range
    char str[128];
    sprintf(str,"Argument above spline interval");
    error->warning(FLERR,str);
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
    // warn if argument below spline range
    char str[128];
    sprintf(str,"Argument below spline interval");
    error->warning(FLERR,str);
  }
  else if(i > coeff_size-1){ 
    i = coeff_size-1;
    // warn if argument above spline range
    char str[128];
    sprintf(str,"Argument above spline interval");
    error->warning(FLERR,str);
  }
  
  double ylo = ystart + i*dy;
  double ybar = (y - ylo)/dy;

  // compute coefficients in y
  
  double a0, a1, a2, a3;
  double p0, p1, p2, p3;
  
  p1 = spline(x,xstart[0],dx[0],coeff[0],coeff_size);
  p2 = spline(x,xstart[1],dx[1],coeff[1],coeff_size);
  
  a0 = p1;

  if(i == 0){
    p3 = spline(x,xstart[2],dx[2],coeff[2],coeff_size);

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
    p3 = spline(x,xstart[2],dx[2],coeff[2],coeff_size);

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
    // warn if argument below spline range
    char str[128];
    sprintf(str,"Argument below spline interval");
    error->warning(FLERR,str);
  }
  else if(i > coeff_size-1){ 
    i = coeff_size-1;
    // warn if argument above spline range
    char str[128];
    sprintf(str,"Argument above spline interval");
    error->warning(FLERR,str);
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
    // warn if argument below spline range
    char str[128];
    sprintf(str,"Argument below spline interval");
    error->warning(FLERR,str);
  }
  else if(i > coeff_size-1){ 
    i = coeff_size-1;
    // warn if argument above spline range
    char str[128];
    sprintf(str,"Argument above spline interval");
    error->warning(FLERR,str);
  }
  
  double ylo = ystart + i*dy;
  double ybar = (y - ylo)/dy;

  // compute coefficients in y
  
  double a0, a1, a2, a3;
  double p0, p1, p2, p3;
  
  p1 = dspline(x,xstart[0],dx[0],coeff[0],coeff_size);
  p2 = dspline(x,xstart[1],dx[1],coeff[1],coeff_size);
  
  a0 = p1;

  if(i == 0){
    p3 = dspline(x,xstart[2],dx[2],coeff[2],coeff_size);

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
    p3 = dspline(x,xstart[2],dx[2],coeff[2],coeff_size);

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
    // warn if argument below spline range
    char str[128];
    sprintf(str,"Argument below spline interval");
    error->warning(FLERR,str);
  }
  else if(i > coeff_size-1){ 
    i = coeff_size-1;
    // warn if argument above spline range
    char str[128];
    sprintf(str,"Argument above spline interval");
    error->warning(FLERR,str);
  }
  
  double ylo = ystart + i*dy;
  double ybar = (y - ylo)/dy;

  // compute coefficients in y
  
  double a0, a1, a2, a3;
  double p0, p1, p2, p3;
  
  p1 = spline(x,xstart[0],dx[0],coeff[0],coeff_size);
  p2 = spline(x,xstart[1],dx[1],coeff[1],coeff_size);
  
  a0 = p1;

  if(i == 0){
    p3 = spline(x,xstart[2],dx[2],coeff[2],coeff_size);

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
    p3 = spline(x,xstart[2],dx[2],coeff[2],coeff_size);

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
      if(i == 0){
        coeff[i][j][0] = data[j][i];
        coeff[i][j][1] = data[j+1][i] - data[j][i];
        coeff[i][j][3] = 0.5*(data[j+2][i] - 2*data[j+1][i] + data[j][i]);
        coeff[i][j][2] = -coeff[i][j][3];
      }
      else if(i == data_size-2){
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

void PairMesoCNT::read_file(char *file, double *data, double *dx, int ninput)
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

  utils::sfgets(FLERR,line,MAXLINE,fp,file,error);
  for(int i = 0; i < ninput; i++){
    if(i > 0) xtemp = x;
    if(NULL == fgets(line,MAXLINE,fp))
      error->one(FLERR,"Premature end of file in pair table");
    if(2 != sscanf(line,"%lg %lg",&x, &data[i])) ++cerror; 
    if(i > 0){
      if(i == 1) *dx = x - xtemp;
      if(*dx != x - xtemp) ++serror;
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
		double *dx, double *dy, int ninput)
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

  utils::sfgets(FLERR,line,MAXLINE,fp,file,error);
  for(int i = 0; i < ninput; i++){
    if(i > 0) ytemp = y;
    for(int j = 0; j < ninput; j++){
      if(j > 0) xtemp = x;
      if(NULL == fgets(line,MAXLINE,fp))
        error->one(FLERR,"Premature end of file in pair table");
      if(3 != sscanf(line,"%lg %lg %lg",&x,&y,&data[j][i])) ++cerror; 
      if(j > 0){
        if(j == 1) dx[i] = x - xtemp;
        if(dx[i] != x - xtemp) ++sxerror;
      }
    }
    if(i > 0){
      if(i == 1) *dy = y - ytemp;
      if(*dy != y - ytemp) ++syerror;
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

double PairMesoCNT::uinf(double h, double alpha, double xi1, double xi2)
{
  double salpha = sin(alpha);
  double salphasq = salpha*salpha;
  if(salphasq < 1.0e-6){
    return (xi2 - xi1) * spline(h,start_uinf,del_uinf,uinf_coeff,pot_points);
  }
  else{
    double omega = 1.0 / (1.0 - comega*salphasq);
    double a = omega * salpha;
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
    double gamma = 1.0 + (gamma_orth - 1.0)*salphasq;

    return gamma * (phi2 - phi1) / a;
  }
}

/* ---------------------------------------------------------------------- */

double PairMesoCNT::usemi(double h, double alpha, 
		double xi1, double xi2, double etaend)
{
  double salpha = sin(alpha);
  double salphasq = salpha*salpha;
  double calpha = cos(alpha);
  double omega = 1.0 / (1.0 - comega*salphasq);
  double theta = 1.0 - ctheta*salphasq; 

  int points = 100;
  double delxi = (xi2 - xi1) / (points -1);
  
  double g, hbar, etabar;
  double sum = 0;

  // first and last term in sum
    
  g = xi1 * omega * salpha;
  hbar = sqrt(h*h + g*g);
  etabar = xi1 * calpha - theta*etaend;
  sum += spline(hbar,etabar,starth_usemi,startxi_usemi,
		  delh_usemi,delxi_usemi,usemi_coeff,pot_points);
  g = xi2 * omega * salpha;
  hbar = sqrt(h*h + g*g);
  etabar = xi2 * calpha - theta*etaend;
  sum += spline(hbar,etabar,starth_usemi,startxi_usemi,
		  delh_usemi,delxi_usemi,usemi_coeff,pot_points);
  sum *= 0.5;

  for(int i = 1; i < points-1; i++){
    double xibar = xi1 + i*delxi;
    g * xibar * omega * salpha;
    hbar = sqrt(h*h + g*g);
    etabar = xibar*calpha - theta*etaend;
    sum += spline(hbar,etabar,starth_usemi,startxi_usemi,
		  delh_usemi,delxi_usemi,usemi_coeff,pot_points);
  }

  double gamma_orth = spline(h,start_gamma,del_gamma,
		    gamma_coeff,gamma_points);
  double gamma = 1.0 + (gamma_orth - 1.0)*salphasq;

  return delxi * gamma * sum;
}

/* ---------------------------------------------------------------------- */

void PairMesoCNT::geom(const double *r1, const double *r2, 
		const double *p1, const double *p2, double *param)
{
  using namespace MathExtra;
  double r[3], p[3], delr[3], l[3], m[3], rbar[3], pbar[3], delrbar[3];
  double psil[3], psim[3], dell_psim[3], delpsil_m[3];
  double delr1[3], delr2[3], delp1[3], delp2[3];
  double ex[3], ey[3];
  double psi, frac, taur, taup;
  double h, alpha, xi1, xi2, eta1, eta2;

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
  frac = 1.0 / (1.0 - psi*psi);
  
  copy3(l,psil);
  scale3(psi,psil);
  copy3(m,psim);
  scale3(psi,psim);

  sub3(l,psim,dell_psim);
  sub3(psil,m,delpsil_m);
  taur = dot3(delr,dell_psim) * frac;
  taup = dot3(delr,delpsil_m) * frac;

  scaleadd3(taur,l,r,rbar);
  scaleadd3(taup,m,p,pbar);
  sub3(pbar,rbar,delrbar);

  h = len3(delrbar);
  
  copy3(delrbar,ex);
  scale3(1/h,ex);
  cross3(l,ex,ey);

  if(dot3(m,ey) < 0) alpha = acos(psi);
  else alpha = MY_2PI - acos(psi);

  sub3(r1,rbar,delr1);
  sub3(r2,rbar,delr2);
  xi1 = dot3(delr1,l);
  xi2 = dot3(delr2,l);

  sub3(p1,pbar,delp1);
  sub3(p2,pbar,delp2);
  eta1 = dot3(delp1,m);
  eta2 = dot3(delp2,m);

  param[0] = h;
  param[1] = alpha;
  param[2] = xi1;
  param[3] = xi2;
  param[4] = eta1;
  param[5] = eta2;
}
