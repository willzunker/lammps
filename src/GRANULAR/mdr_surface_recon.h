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
   MDR deformed-particle surface reconstruction.

   Faithful C++ port of the MATLAB reconstruction (capmaker.m / spheremaker.m
   / plot_particles3.m) used with the MDR contact model.  Each contact carves
   a flattened "cap" (surface of revolution about the contact normal); the
   deformed particle is the per-direction minimum over all caps and the
   residual sphere of apparent radius R.

   This header has no LAMMPS dependencies so it can be unit-tested standalone.

   Contributing authors:
   Dalil Ashong (UC Berkeley), William Zunker (MIT), Ken Kamrin (UC Berkeley)
------------------------------------------------------------------------- */

#ifndef LMP_MDR_SURFACE_RECON_H
#define LMP_MDR_SURFACE_RECON_H

#include <algorithm>
#include <cmath>
#include <vector>

namespace LAMMPS_NS {
namespace Granular_MDR_NS {

  namespace recon_detail {
    constexpr double RPI = 3.14159265358979323846;
    constexpr double RHALFPI = 1.57079632679489661923;
    constexpr double RTWOPI = 6.28318530717958647692;
  }    // namespace recon_detail

  // Complete elliptic integral of the second kind, MATLAB ellipticE(m)
  // convention (parameter m = k^2), valid for m in [0,1].  AGM method,
  // machine precision.
  static inline double mdr_ellipticE(double m)
  {
    using namespace recon_detail;
    if (m <= 0.0) return RHALFPI;
    if (m >= 1.0) return 1.0;
    double a = 1.0, b = std::sqrt(1.0 - m), c = std::sqrt(m);
    double s = 0.5 * m;       // n = 0 term: 2^{-1} c0^2
    double two_pow = 1.0;     // n = 1 term weight 2^{0}, doubling thereafter
    for (int i = 0; i < 80; i++) {
      const double an = 0.5 * (a + b);
      const double bn = std::sqrt(a * b);
      const double cn = 0.5 * (a - b);
      a = an; b = bn; c = cn;
      s += two_pow * c * c;
      two_pow *= 2.0;
      if (std::fabs(c) < 1e-17) break;
    }
    const double K = RPI / (2.0 * a);
    return K * (1.0 - s);
  }

  // One contact's cap as a meridian profile in the contact-local frame
  // (polar angle alpha measured from the contact normal n).  Sampled so that
  // alpha is strictly increasing; r is the distance from the particle centre.
  // For alpha > alpha_max the cap imposes no constraint (residual sphere).
  struct MDRCapProfile {
    std::vector<double> alpha;    // increasing, radians, starts at 0
    std::vector<double> r;        // radius from centre at each alpha
    double alpha_max = 0.0;       // angular extent of the cap footprint

    // radius this cap allows along a ray at polar angle a from n;
    // returns a huge value outside the footprint (no constraint)
    double radius_at(double a) const
    {
      if (a >= alpha_max) return 1e30;
      // binary search in increasing alpha
      const auto it = std::upper_bound(alpha.begin(), alpha.end(), a);
      size_t hi = it - alpha.begin();
      if (hi == 0) return r.front();
      if (hi >= alpha.size()) return r.back();
      const size_t lo = hi - 1;
      const double t = (a - alpha[lo]) / (alpha[hi] - alpha[lo]);
      return r[lo] + t * (r[hi] - r[lo]);
    }
  };

  // append a (rho,z) meridian sample, converting to (alpha,r); keeps alpha
  // strictly increasing (drops non-monotone duplicates from numerical noise)
  static inline void recon_push(MDRCapProfile &p, double rho, double z)
  {
    const double r = std::sqrt(rho * rho + z * z);
    const double a = std::atan2(rho, z);
    if (!p.alpha.empty() && a <= p.alpha.back() + 1e-14) return;
    p.alpha.push_back(a);
    p.r.push_back(r);
  }

  // sanity check on a built meridian; rejects degenerate/over-carved/blown-up
  // profiles so a clean plane-clip fallback can be substituted. rfloor is the
  // contact-plane distance R-delta: a physical cap never carves nearer the centre
  // than the contact plane.
  static inline bool recon_profile_valid(const MDRCapProfile &p, double R, double rfloor)
  {
    if (p.alpha.size() < 4) return false;
    if (p.alpha.back() > 2.7) return false;    // cap footprint > ~155 deg is unphysical
    for (double rv : p.r)
      if (!(rv >= rfloor && rv < 1.7 * R)) return false;
    return true;
  }

  // plane-clip cap: sphere of radius R cut by a plane a distance zflat from the
  // centre (toward the contact). Robustness fallback for degenerate transitions.
  static inline MDRCapProfile mdr_plane_clip_profile(double R, double zflat, int nseg = 200)
  {
    MDRCapProfile p;
    if (zflat < R && zflat > 0.0) {
      const double aclip = std::sqrt(std::max(0.0, R * R - zflat * zflat));
      for (int i = 0; i <= nseg; i++) recon_push(p, aclip * i / nseg, zflat);
    }
    p.alpha_max = p.alpha.empty() ? 0.0 : p.alpha.back();
    return p;
  }

  // Lift a cap meridian outward along the contact normal by dz (toward the
  // surface), used to apply the elastic spring-back to a persisted plastic
  // imprint after the contact separates. Points lifted past the sphere are
  // clipped later by the min-over-caps against R.
  static inline MDRCapProfile mdr_lift_profile(const MDRCapProfile &in, double dz)
  {
    MDRCapProfile p;
    for (size_t i = 0; i < in.alpha.size(); i++) {
      const double rho = in.r[i] * std::sin(in.alpha[i]);
      const double z = in.r[i] * std::cos(in.alpha[i]) + dz;
      recon_push(p, rho, z);
    }
    p.alpha_max = p.alpha.empty() ? 0.0 : p.alpha.back();
    return p;
  }

  // Build the cap meridian profile for a single contact side.
  //   yflag : 0 elastic, otherwise plastic
  //   R     : apparent radius
  //   delta : (partitioned) apparent overlap on this side
  //   deltae: cap "deltae" field (deltamax_MDR - deltap)   [capmaker deltae1]
  //   A,B   : elliptical indenter geometry
  //   a_na  : non-adhesive contact radius   [capmaker a1]
  //   deltamax : maximum MDR apparent overlap on this side [capmaker deltamax1]
  //   deltaR   : plastic displacement correction
  //   amax     : maximum experienced contact radius
  //   nseg     : meridian samples per segment
  static inline MDRCapProfile mdr_cap_profile(double yflag, double R, double delta, double deltae,
                                              double A, double B, double a_na, double deltamax,
                                              double deltaR, double amax, int nseg = 200)
  {
    using namespace recon_detail;
    MDRCapProfile p;
    const double zflat = R - delta;    // contact-face plane distance from centre

    if (yflag == 0.0) {
      // ---- elastic: flat contact face then sphere ----
      const double thetafin = std::asin(std::max(-1.0, std::min(1.0, (R - delta - deltae) / R)));
      const double xfin = R * std::cos(thetafin);    // rho where cap rejoins sphere
      const double a1 = std::min(a_na, xfin);

      for (int i = 0; i <= nseg; i++) {
        const double rho = a1 * i / nseg;
        recon_push(p, rho, zflat);
      }
      for (int i = 1; i <= nseg; i++) {
        const double rho = a1 + (xfin - a1) * i / nseg;
        const double z = std::sqrt(std::max(0.0, R * R - rho * rho));
        recon_push(p, rho, z);
      }
    } else {
      // ---- plastic: flat face then scaled transition (elliptic + parabola) ----
      const double deltaTR = (delta - deltamax + 0.5 * A + deltaR) / (1.0 + (0.5 * A) / deltaR);
      const double thetafin =
          std::asin(std::max(-1.0, std::min(1.0, (R - deltamax - deltae) / R)));
      const double yfin = R * std::sin(thetafin);

      const double h2 = R - deltamax - (0.5 * A + deltaR);
      const double x2 = std::sqrt(std::max(0.0, R * R - h2 * h2));
      const double x1 = amax;

      // parabola a x^2 + b x + c through (x1, yell_end), (x2, h2), slope 0 at x1
      const double mend = 4.0 * amax * amax / (B * B) - 0.0005;
      const double yell_end = (R - delta - deltaTR) - (0.5 / RPI) * A * (RPI - 2.0 * mdr_ellipticE(mend));
      const double det = (x1 * x1) * (x2 - 1.0) - x1 * (x2 * x2 - 1.0) + (x2 * x2 - x2);
      // solve 3x3 directly: rows [x1^2 x1 1; x2^2 x2 1; 2x1 1 0]
      double pa, pb, pc;
      {
        const double a11 = x1 * x1, a12 = x1, a13 = 1.0;
        const double a21 = x2 * x2, a22 = x2, a23 = 1.0;
        const double a31 = 2.0 * x1, a32 = 1.0, a33 = 0.0;
        const double v1 = yell_end, v2 = h2, v3 = 0.0;
        const double D = a11 * (a22 * a33 - a23 * a32) - a12 * (a21 * a33 - a23 * a31) +
            a13 * (a21 * a32 - a22 * a31);
        pa = (v1 * (a22 * a33 - a23 * a32) - a12 * (v2 * a33 - a23 * v3) +
              a13 * (v2 * a32 - a22 * v3)) / D;
        pb = (a11 * (v2 * a33 - a23 * v3) - v1 * (a21 * a33 - a23 * a31) +
              a13 * (a21 * v3 - v2 * a31)) / D;
        pc = (a11 * (a22 * v3 - v2 * a32) - a12 * (a21 * v3 - v2 * a31) +
              v1 * (a21 * a32 - a22 * a31)) / D;
        (void) det;
      }

      auto fullprof = [&](double rho) -> double {
        if (rho < amax) {
          const double ms = 4.0 * rho * rho / (B * B);
          return (R - delta - deltaTR) - (0.5 / RPI) * A * (RPI - 2.0 * mdr_ellipticE(ms));
        }
        return pa * rho * rho + pb * rho + pc;
      };

      const double fp_a = fullprof(a_na);
      const double scaler = ((R - delta) - yfin) / (fp_a - yfin);

      // capmaker's scaled transition is only well-posed for a loaded contact with
      // a finite contact radius. For (near-)unloaded contacts (a_na -> 0) or an
      // inverted stitch (scaler <= 0) it degenerates into spurious egg / saw-tooth
      // shapes, so detect that and fall back to a clean plane clip below.
      const bool ok =
          (a_na > 1e-9 * R) && (x2 > a_na) && std::isfinite(scaler) && (scaler > 0.0);
      if (ok) {
        for (int i = 0; i <= nseg; i++) recon_push(p, a_na * i / nseg, zflat);
        for (int i = 1; i <= nseg; i++) {
          const double rho = a_na + (x2 - a_na) * i / nseg;
          recon_push(p, rho, scaler * (fullprof(rho) - yfin) + yfin);
        }
      }
    }

    // Universal robustness guard (both branches): extreme or unloaded contacts in
    // dense compaction can invert capmaker's scaled transition into spurious
    // egg/saw-tooth shapes that carve below the contact plane. Reject any such
    // meridian and substitute a clean plane clip at the contact plane z = R-delta.
    if (!p.alpha.empty() && !recon_profile_valid(p, R, 0.93 * zflat)) {
      p.alpha.clear();
      p.r.clear();
    }
    if (p.alpha.empty()) return mdr_plane_clip_profile(R, zflat, nseg);

    p.alpha_max = p.alpha.back();
    return p;
  }

  // one contact contributing to a particle: outward normal + cap profile
  struct MDRCap {
    double n[3];
    MDRCapProfile profile;
  };

  // deformed-surface radius along unit ray d (from particle centre), as the
  // minimum over all caps and the residual sphere of apparent radius R
  static inline double mdr_surface_radius(const double d[3], double R,
                                          const std::vector<MDRCap> &caps)
  {
    double rmin = R;
    for (const auto &c : caps) {
      double dot = d[0] * c.n[0] + d[1] * c.n[1] + d[2] * c.n[2];
      if (dot > 1.0) dot = 1.0;
      if (dot < -1.0) dot = -1.0;
      const double a = std::acos(dot);
      if (a < c.profile.alpha_max) {
        const double rc = c.profile.radius_at(a);
        if (rc < rmin) rmin = rc;
      }
    }
    return rmin;
  }

  // triangulated deformed-particle surface as a UV sphere (shared poles).
  // ntheta = azimuthal divisions, nphi = polar divisions (>= 2).
  struct MDRMesh {
    std::vector<double> xyz;    // 3 per vertex (lab frame, includes centre)
    std::vector<int> tris;     // 3 vertex indices per triangle
  };

  static inline void mdr_ray(double phi, double theta, double d[3])
  {
    const double sp = std::sin(phi);
    d[0] = sp * std::cos(theta);
    d[1] = sp * std::sin(theta);
    d[2] = std::cos(phi);
  }

  static inline MDRMesh mdr_build_particle_mesh(const double center[3], double R,
                                                const std::vector<MDRCap> &caps, int ntheta,
                                                int nphi)
  {
    using namespace recon_detail;
    MDRMesh m;
    if (ntheta < 3) ntheta = 3;
    if (nphi < 2) nphi = 2;

    auto vertex = [&](double phi, double theta) {
      double d[3];
      mdr_ray(phi, theta, d);
      const double r = mdr_surface_radius(d, R, caps);
      m.xyz.push_back(center[0] + r * d[0]);
      m.xyz.push_back(center[1] + r * d[1]);
      m.xyz.push_back(center[2] + r * d[2]);
    };

    // north pole, (nphi-1) interior rings, south pole
    const int nring = nphi - 1;
    const int north = 0;
    vertex(0.0, 0.0);
    const int ring0 = 1;
    for (int i = 0; i < nring; i++) {
      const double phi = RPI * (i + 1) / nphi;
      for (int j = 0; j < ntheta; j++) vertex(phi, RTWOPI * j / ntheta);
    }
    const int south = 1 + nring * ntheta;
    vertex(RPI, 0.0);

    auto ringidx = [&](int i, int j) { return ring0 + i * ntheta + (j % ntheta); };

    // north fan
    for (int j = 0; j < ntheta; j++) {
      m.tris.push_back(north);
      m.tris.push_back(ringidx(0, j));
      m.tris.push_back(ringidx(0, j + 1));
    }
    // bands between rings
    for (int i = 0; i < nring - 1; i++)
      for (int j = 0; j < ntheta; j++) {
        const int a = ringidx(i, j), b = ringidx(i, j + 1);
        const int c = ringidx(i + 1, j), d = ringidx(i + 1, j + 1);
        m.tris.push_back(a); m.tris.push_back(c); m.tris.push_back(b);
        m.tris.push_back(b); m.tris.push_back(c); m.tris.push_back(d);
      }
    // south fan
    for (int j = 0; j < ntheta; j++) {
      m.tris.push_back(south);
      m.tris.push_back(ringidx(nring - 1, j + 1));
      m.tris.push_back(ringidx(nring - 1, j));
    }
    return m;
  }

}    // namespace Granular_MDR_NS
}    // namespace LAMMPS_NS

#endif
