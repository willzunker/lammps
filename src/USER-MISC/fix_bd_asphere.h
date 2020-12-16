/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(bd/asphere,FixBdAsphere)

#else

#ifndef LMP_FIX_BD_ASPHERE_H
#define LMP_FIX_BD_ASPHERE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixBdAsphere : public Fix {
 public:
  FixBdAsphere(class LAMMPS *, int, char **);
  virtual ~FixBdAsphere();
  void init();
  void initial_integrate(int);
  void setup(int);
  void post_force(int);
  int setmask();
  void reset_dt();
  void update_x_and_omega(double *, double *, double *,
			  double *, double *, int );


 private:
  int seed;               // RNG seed
  int dipole_flag;        // set if dipole is used
  double dt, sqrtdt;      // time step interval and its sqrt


  double gamma_t,gamma_r;  // translational and rotational damping params
  double diff_t,diff_r;    // translational and rotational diffusion coeffs

  double g1,g2, g3, g4;    // prefactors in time stepping
  int noise_flag;          // 0/1 for noise off/on
  int gaussian_noise_flag; // 0/1 for uniform/gaussian noise

  class AtomVecEllipsoid *avec;
  
protected:
  class RanMars *random;
  typedef double (RanMars::*rng_member)();
  rng_member rng_func;    // placeholder for RNG function

};

}
#endif
#endif

/* ERROR/WARNING messages:

E: Illegal fix bd/asphere command.

Wrong number/type of input arguments.

E: Compute bd/asphere requires atom style sphere

Self-explanatory.

E: Compute bd/asphere requires atom style ellipsoid

Self-explanatory.

E: Compute bd/asphere dipole requires atom attribute mu

Self-explanatory.

E: Fix bd/asphere translational viscous drag coefficient must be > 0.

Self-explanatory.

E: Fix bd/asphere rotational viscous drag coefficient must be > 0.

Self-explanatory.

E: Fix bd/asphere translational diffusion coefficient must be > 0.

Self-explanatory.

E: Fix bd/asphere rotational diffusion coefficient must be > 0.

Self-explanatory.

E: Fix bd/asphere seed must be > 0.

*/
