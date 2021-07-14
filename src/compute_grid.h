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

#ifndef LMP_COMPUTE_GRID_H
#define LMP_COMPUTE_GRID_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeGrid : public Compute {
 public:

  ComputeGrid(class LAMMPS *, int, char **);
  virtual ~ComputeGrid();
  void init();
  void setup();
  virtual void compute_array() = 0;

  double memory_usage();

 protected:
  int nx, ny, nz;                      // global grid dimensions
  int nxlo, nxhi, nylo, nyhi, nzlo, nzhi; // local grid bounds, inclusive
  int ngrid;                           // number of global grid points
  int ngridlocal;                      // number of local grid points
  int nvalues;                         // number of values per grid point
  double **grid;                       // global grid
  double **gridall;                    // global grid summed over procs
  double ****gridlocal;                // local grid
  int triclinic;                       // triclinic flag
  double *boxlo, *prd;                 // box info (units real/ortho or reduced/tri)
  double *sublo, *subhi;               // subdomain info (units real/ortho or reduced/tri)
  double delxinv,delyinv,delzinv;      // inverse grid spacing
  double delx,dely,delz;               // grid spacing
  int nargbase;                        // number of base class args 
  double cutmax;                       // largest cutoff distance
  int size_array_cols_base;            // number of columns used for coords, etc.
  int *local_flags;                    // local flag for each grid point
  int gridlocal_allocated;             // shows if gridlocal allocated

  void allocate();
  void grid2x(int, double*);           // convert grid point to coord
  void grid2ix(int, int&, int&, int&); // convert grid point to ix, iy, iz
  void assign_coords();                // assign coords for grid
  void assign_local_flags();           // set local flag for each grid point
  int check_local(int);                // check if grid point is local
  void set_grid_global();              // set global grid
  void set_grid_local();               // set bounds for local grid
 private:
};

}

#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
