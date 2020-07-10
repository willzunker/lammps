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

#ifndef LMP_MLIAP_MODEL_LINEAR_H
#define LMP_MLIAP_MODEL_LINEAR_H

#include "mliap_model.h"

namespace LAMMPS_NS {

class MLIAPModelLinear : public MLIAPModel {
public:
  MLIAPModelLinear(LAMMPS*, char*);
  MLIAPModelLinear(LAMMPS*, int, int);
  ~MLIAPModelLinear();
  virtual void gradient(int, int*, int*, double**, double**, class PairMLIAP*, int);
  virtual void param_gradient(int, int*, int*, double**, int**, int**, double**, double*);
  virtual int get_gamma_nnz();
  virtual void compute_force_gradients(double**, int, int*, int*, int*, int*, 
                                       int*, double***, int, int, double**, double*);
protected:
};

}

#endif

