/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(wall/gran/kk,FixWallGranKokkos<LMPDeviceType>)
FixStyle(wall/gran/kk/device,FixWallGranKokkos<LMPDeviceType>)
FixStyle(wall/gran/kk/host,FixWallGranKokkos<LMPHostType>)

#else

#ifndef LMP_FIX_WALL_GRAN_KOKKOS_H
#define LMP_FIX_WALL_GRAN_KOKKOS_H

#include "fix_wall_gran.h"
#include "kokkos_type.h"
#include "kokkos_base.h"

namespace LAMMPS_NS {

template<class DeviceType>
class FixWallGranKokkos : public FixWallGran, public KokkosBase {
 public:
  FixWallGranKokkos(class LAMMPS *, int, char **);
  ~FixWallGranKokkos();
  void init();
  void post_force(int);
  void grow_arrays(int);
  void copy_arrays(int, int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  int pack_exchange_kokkos(const int &nsend,DAT::tdual_xfloat_2d &buf,
                           DAT::tdual_int_1d k_sendlist,
                           DAT::tdual_int_1d k_copylist,
                           ExecutionSpace space, int dim,
                           X_FLOAT lo, X_FLOAT hi);
  void unpack_exchange_kokkos(DAT::tdual_xfloat_2d &k_buf,
                              DAT::tdual_int_1d &indices,int nrecv,
                              int nlocal,int dim,X_FLOAT lo,X_FLOAT hi,
                              ExecutionSpace space);

  template <int WallStyle>
  KOKKOS_INLINE_FUNCTION
  void hooke_history_item(const int &i) const;

 protected:
  X_FLOAT wlo;
  X_FLOAT whi;
  V_FLOAT vwall[3];

  typedef ArrayTypes<DeviceType> AT;
  typename AT::t_x_array x;
  typename AT::t_v_array v;
  typename AT::t_v_array omega_;
  typename AT::t_f_array f;
  typename AT::t_f_array torque;
  typename AT::t_int_1d mask;
  typename AT::t_float_1d rmass;
  typename AT::t_float_1d radius_;
  typename AT::tdual_float_2d k_history_one;
  typename AT::t_float_2d d_history_one;
};

template <class DeviceType, int WallStyle>
struct FixWallGranKokkosHookeHistoryFunctor {
  FixWallGranKokkos<DeviceType> c;
  FixWallGranKokkosHookeHistoryFunctor(FixWallGranKokkos<DeviceType> *c_ptr): c(*c_ptr) {}
  KOKKOS_INLINE_FUNCTION
  void operator()(const int &i) const {
    c.template hooke_history_item<WallStyle>(i);
  }
};
}

#endif
#endif
