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

/* ----------------------------------------------------------------------
   Shared MDR contact geometry.

   The deformed-surface ("cap") parameters below are computed identically
   for the force evaluation (GranSubModNormalMDR::calculate_forces) and for
   surface reconstruction (dump mdr/surface), so both consume a single
   source of truth.  The expressions are copied verbatim from the force
   routine; keep them in lock-step.

   Contributing authors:
   Dalil Ashong (UC Berkeley), William Zunker (MIT), Ken Kamrin (UC Berkeley)
------------------------------------------------------------------------- */

#ifndef LMP_MDR_RECONSTRUCT_H
#define LMP_MDR_RECONSTRUCT_H

#include "math_const.h"

#include <cmath>

namespace LAMMPS_NS {
namespace Granular_MDR_NS {

  // average pressure along the yield surface
  static inline double mdr_pressure_yield(double Y, double deltamax_MDR, double R)
  {
    return Y * (1.75 * exp(-4.4 * deltamax_MDR / R) + 1.0);
  }

  // deformed-surface (cap) geometry for one side of one contact
  struct MDRCapGeometry {
    double A;          // height of elliptical indenter
    double Ainv;       // 1/A
    double B;          // width of elliptical indenter
    double deltae1D;   // transformed elastic displacement
    double deltaR;     // displacement correction (0 in the elastic regime)
    double amax;       // maximum experienced contact radius
    double a_na;       // non-adhesive contact radius
    double deltap;     // rigid-flat offset (plastic regime; 0 otherwise)
  };

  // Reproduces the cap geometry block of GranSubModNormalMDR::calculate_forces.
  // yflag: 0 elastic, otherwise plastic.  pY from mdr_pressure_yield().
  static inline MDRCapGeometry mdr_cap_geometry(double yflag, double R, double delta_MDR,
                                                double deltamax_MDR, double cA, double pY,
                                                double Eeff, double Eeffinv, double G, double poiss)
  {
    constexpr double PIINV = 0.318309886183790691216;    // 1/PI
    using MathConst::MY_2PI;

    MDRCapGeometry c;
    c.deltaR = 0.0;
    c.deltap = 0.0;

    double A, Ainv, B, deltae1D, deltaR, amax, amaxsq;

    if (yflag == 0.0) {
      // elastic contact
      A = 4.0 * R;
      Ainv = 1.0 / A;
      B = 2.0 * R;
      deltae1D = delta_MDR;
      amax = sqrt(deltamax_MDR * R);
    } else {
      // plastic contact
      amax = sqrt(2.0 * deltamax_MDR * R - pow(deltamax_MDR, 2) + cA * PIINV);
      amaxsq = amax * amax;
      A = 4.0 * pY * Eeffinv * amax;
      Ainv = 1.0 / A;
      B = 2.0 * amax;

      // maximum transformed elastic displacement
      const double deltae1Dmax = A * 0.5;

      // force caused by full submersion of elliptical indenter to depth of A/2
      double Fmax = Eeff * (A * B * 0.25) * acos(1 - 2 * deltae1Dmax * Ainv);
      Fmax -= (2 - 4 * deltae1Dmax * Ainv) * sqrt(deltae1Dmax * Ainv - pow(deltae1Dmax * Ainv, 2));

      // depth of particle center
      const double zR = R - (deltamax_MDR - deltae1Dmax);

      deltaR = 2 * amaxsq * (-1 + poiss) - (-1 + 2 * poiss) * zR * (-zR + sqrt(amaxsq + pow(zR, 2)));
      deltaR *= Fmax / (MY_2PI * amaxsq * G * sqrt(amaxsq + pow(zR, 2)));

      // transformed elastic displacement
      deltae1D = (delta_MDR - deltamax_MDR + deltae1Dmax + deltaR) / (1 + deltaR / deltae1Dmax);

      c.deltaR = deltaR;

      // added for rigid flat placement
      c.deltap = deltamax_MDR - (deltae1Dmax + deltaR);
    }

    double a_na;
    (deltae1D >= 0.0) ? a_na = B * sqrt(A - deltae1D) * sqrt(deltae1D) * Ainv : a_na = 0.0;

    c.A = A;
    c.Ainv = Ainv;
    c.B = B;
    c.deltae1D = deltae1D;
    c.amax = amax;
    c.a_na = a_na;
    return c;
  }

}    // namespace Granular_MDR_NS
}    // namespace LAMMPS_NS

#endif
