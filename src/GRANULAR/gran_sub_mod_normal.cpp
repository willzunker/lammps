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

#include "gran_sub_mod_normal.h"
#include "error.h"
#include "granular_model.h"
#include "math_const.h"
#include "atom.h"
#include "csv_writer.h"
#include "update.h"

#include <cmath>
#include <iostream>
#include <iomanip> 
#include <sstream>

using namespace LAMMPS_NS;
using namespace Granular_NS;

using MathConst::MY_2PI;
using MathConst::MY_PI;

static constexpr double PI27SQ = 266.47931882941264802866;      // 27*PI**2
static constexpr double THREEROOT3 = 5.19615242270663202362;    // 3*sqrt(3)
static constexpr double SIXROOT6 = 14.69693845669906728801;     // 6*sqrt(6)
static constexpr double INVROOT6 = 0.40824829046386307274;      // 1/sqrt(6)
static constexpr double FOURTHIRDS = (4.0 / 3.0);               // 4/3
static constexpr double JKRPREFIX = 1.2277228507842888;         // cbrt(3*PI**2/16)

/* ----------------------------------------------------------------------
   Default normal model
------------------------------------------------------------------------- */

GranSubModNormal::GranSubModNormal(GranularModel *gm, LAMMPS *lmp) : GranSubMod(gm, lmp)
{
  material_properties = 0;
  cohesive_flag = 0;
}

/* ---------------------------------------------------------------------- */

bool GranSubModNormal::touch()
{
  bool touchflag = (gm->rsq < gm->radsum * gm->radsum);
  return touchflag;
}

/* ---------------------------------------------------------------------- */

double GranSubModNormal::pulloff_distance(double /*radi*/, double /*radj*/)
{
  // called outside of compute(), do not assume correct geometry defined in contact
  return 0.0;
}

/* ---------------------------------------------------------------------- */

double GranSubModNormal::calculate_contact_radius()
{
  return sqrt(gm->dR);
}

/* ---------------------------------------------------------------------- */

void GranSubModNormal::set_fncrit()
{
  Fncrit = fabs(gm->Fntot);
}

/* ----------------------------------------------------------------------
   No model
------------------------------------------------------------------------- */

GranSubModNormalNone::GranSubModNormalNone(GranularModel *gm, LAMMPS *lmp) :
    GranSubModNormal(gm, lmp)
{
}

/* ---------------------------------------------------------------------- */

double GranSubModNormalNone::calculate_forces()
{
  return 0.0;
}

/* ----------------------------------------------------------------------
   Hookean normal force
------------------------------------------------------------------------- */

GranSubModNormalHooke::GranSubModNormalHooke(GranularModel *gm, LAMMPS *lmp) :
    GranSubModNormal(gm, lmp)
{
  num_coeffs = 2;
}

/* ---------------------------------------------------------------------- */

void GranSubModNormalHooke::coeffs_to_local()
{
  k = coeffs[0];
  damp = coeffs[1];

  if (k < 0.0 || damp < 0.0) error->all(FLERR, "Illegal Hooke normal model");
}

/* ---------------------------------------------------------------------- */

double GranSubModNormalHooke::calculate_forces()
{
  return k * gm->delta;
}

/* ----------------------------------------------------------------------
   Hertzian normal force
------------------------------------------------------------------------- */

GranSubModNormalHertz::GranSubModNormalHertz(GranularModel *gm, LAMMPS *lmp) :
    GranSubModNormal(gm, lmp)
{
  num_coeffs = 2;
  contact_radius_flag = 1;
}

/* ---------------------------------------------------------------------- */

void GranSubModNormalHertz::coeffs_to_local()
{
  k = coeffs[0];
  damp = coeffs[1];

  if (k < 0.0 || damp < 0.0) error->all(FLERR, "Illegal Hertz normal model");
}

/* ---------------------------------------------------------------------- */

double GranSubModNormalHertz::calculate_forces()
{
  return k * gm->contact_radius * gm->delta;
}

/* ----------------------------------------------------------------------
   Hertzian normal force with material properties
------------------------------------------------------------------------- */

GranSubModNormalHertzMaterial::GranSubModNormalHertzMaterial(GranularModel *gm, LAMMPS *lmp) :
    GranSubModNormalHertz(gm, lmp)
{
  material_properties = 1;
  num_coeffs = 3;
  contact_radius_flag = 1;
  mixed_coefficients = 0;
}

/* ---------------------------------------------------------------------- */

void GranSubModNormalHertzMaterial::coeffs_to_local()
{
  Emod = coeffs[0];
  damp = coeffs[1];
  poiss = coeffs[2];
  if (!mixed_coefficients) {
    if (gm->contact_type == PAIR) {
      k = FOURTHIRDS * mix_stiffnessE(Emod, Emod, poiss, poiss);
    } else {
      k = FOURTHIRDS * mix_stiffnessE_wall(Emod, poiss);
    }
  }

  if (Emod < 0.0 || damp < 0.0) error->all(FLERR, "Illegal Hertz material normal model");
}

/* ---------------------------------------------------------------------- */

void GranSubModNormalHertzMaterial::mix_coeffs(double *icoeffs, double *jcoeffs)
{
  coeffs[0] = mix_stiffnessE(icoeffs[0], jcoeffs[0], icoeffs[2], jcoeffs[2]);
  coeffs[1] = mix_geom(icoeffs[1], jcoeffs[1]);
  coeffs[2] = mix_geom(icoeffs[2], jcoeffs[2]);

  k = FOURTHIRDS * coeffs[0];
  mixed_coefficients = 1;

  coeffs_to_local();
}

/* ----------------------------------------------------------------------
   DMT normal force
------------------------------------------------------------------------- */

GranSubModNormalDMT::GranSubModNormalDMT(GranularModel *gm, LAMMPS *lmp) : GranSubModNormal(gm, lmp)
{
  material_properties = 1;
  cohesive_flag = 1;
  num_coeffs = 4;
  contact_radius_flag = 1;
  mixed_coefficients = 0;
}

/* ---------------------------------------------------------------------- */

void GranSubModNormalDMT::coeffs_to_local()
{
  Emod = coeffs[0];
  damp = coeffs[1];
  poiss = coeffs[2];
  cohesion = coeffs[3];

  if (!mixed_coefficients) {
    if (gm->contact_type == PAIR) {
      k = FOURTHIRDS * mix_stiffnessE(Emod, Emod, poiss, poiss);
    } else {
      k = FOURTHIRDS * mix_stiffnessE_wall(Emod, poiss);
    }
  }

  if (Emod < 0.0 || damp < 0.0) error->all(FLERR, "Illegal DMT normal model");
}

/* ---------------------------------------------------------------------- */

void GranSubModNormalDMT::mix_coeffs(double *icoeffs, double *jcoeffs)
{
  coeffs[0] = mix_stiffnessE(icoeffs[0], jcoeffs[0], icoeffs[2], jcoeffs[2]);
  coeffs[1] = mix_geom(icoeffs[1], jcoeffs[1]);
  coeffs[2] = mix_geom(icoeffs[2], jcoeffs[2]);
  coeffs[3] = mix_geom(icoeffs[3], jcoeffs[3]);

  k = FOURTHIRDS * coeffs[0];
  mixed_coefficients = 1;

  coeffs_to_local();
}

/* ---------------------------------------------------------------------- */

double GranSubModNormalDMT::calculate_forces()
{
  Fne = k * gm->contact_radius * gm->delta;
  F_pulloff = 4.0 * MY_PI * cohesion * gm->Reff;
  Fne -= F_pulloff;
  return Fne;
}

/* ---------------------------------------------------------------------- */

void GranSubModNormalDMT::set_fncrit()
{
  Fncrit = fabs(Fne + 2.0 * F_pulloff);
}

/* ----------------------------------------------------------------------
   JKR normal force
------------------------------------------------------------------------- */

GranSubModNormalJKR::GranSubModNormalJKR(GranularModel *gm, LAMMPS *lmp) : GranSubModNormal(gm, lmp)
{
  material_properties = 1;
  cohesive_flag = 1;
  beyond_contact = 1;
  num_coeffs = 4;
  contact_radius_flag = 1;
  mixed_coefficients = 0;
}

/* ---------------------------------------------------------------------- */

void GranSubModNormalJKR::coeffs_to_local()
{
  Emod = coeffs[0];
  damp = coeffs[1];
  poiss = coeffs[2];
  cohesion = coeffs[3];

  if (!mixed_coefficients) {
    if (gm->contact_type == PAIR) {
      Emix = mix_stiffnessE(Emod, Emod, poiss, poiss);
    } else {
      Emix = mix_stiffnessE_wall(Emod, poiss);
    }
  }

  k = FOURTHIRDS * Emix;

  if (Emod < 0.0 || damp < 0.0) error->all(FLERR, "Illegal JKR normal model");
}

/* ---------------------------------------------------------------------- */

void GranSubModNormalJKR::mix_coeffs(double *icoeffs, double *jcoeffs)
{
  coeffs[0] = mix_stiffnessE(icoeffs[0], jcoeffs[0], icoeffs[2], jcoeffs[2]);
  coeffs[1] = mix_geom(icoeffs[1], jcoeffs[1]);
  coeffs[2] = mix_geom(icoeffs[2], jcoeffs[2]);
  coeffs[3] = mix_geom(icoeffs[3], jcoeffs[3]);

  Emix = coeffs[0];
  mixed_coefficients = 1;

  coeffs_to_local();
}

/* ---------------------------------------------------------------------- */

bool GranSubModNormalJKR::touch()
{
  double delta_pulloff, dist_pulloff;
  bool touchflag;

  if (gm->touch) {
    // delta_pulloff defined as positive so center-to-center separation is > radsum
    delta_pulloff = JKRPREFIX * cbrt(gm->Reff * cohesion * cohesion / (Emix * Emix));
    dist_pulloff = gm->radsum + delta_pulloff;
    touchflag = gm->rsq < (dist_pulloff * dist_pulloff);
  } else {
    touchflag = gm->rsq < (gm->radsum * gm->radsum);
  }

  return touchflag;
}

/* ----------------------------------------------------------------------
  called outside of compute(), do not assume geometry defined in contact
------------------------------------------------------------------------- */

double GranSubModNormalJKR::pulloff_distance(double radi, double radj)
{
  double Reff_tmp;

  Reff_tmp = radi * radj / (radi + radj);    // May not be defined
  if (Reff_tmp <= 0) return 0;
  // Defined as positive so center-to-center separation is > radsum
  return JKRPREFIX * cbrt(Reff_tmp * cohesion * cohesion / (Emix * Emix));
}

/* ---------------------------------------------------------------------- */

double GranSubModNormalJKR::calculate_contact_radius()
{
  double R2, dR2, t0, t1, t2, t3, t4, t5, t6;
  double sqrt1, sqrt2, sqrt3;

  R2 = gm->Reff * gm->Reff;
  dR2 = gm->dR * gm->dR;
  t0 = cohesion * cohesion * R2 * R2 * Emix;
  t1 = PI27SQ * t0;
  t2 = 8.0 * gm->dR * dR2 * Emix * Emix * Emix;
  t3 = 4.0 * dR2 * Emix;

  // in case sqrt(0) < 0 due to precision issues
  sqrt1 = MAX(0, t0 * (t1 + 2.0 * t2));
  t4 = cbrt(t1 + t2 + THREEROOT3 * MY_PI * sqrt(sqrt1));
  t5 = t3 / t4 + t4 / Emix;
  sqrt2 = MAX(0, 2.0 * gm->dR + t5);
  t6 = sqrt(sqrt2);
  sqrt3 = MAX(0, 4.0 * gm->dR - t5 + SIXROOT6 * cohesion * MY_PI * R2 / (Emix * t6));

  return INVROOT6 * (t6 + sqrt(sqrt3));
}

/* ---------------------------------------------------------------------- */

double GranSubModNormalJKR::calculate_forces()
{
  double a2;
  a2 = gm->contact_radius * gm->contact_radius;
  Fne = k * gm->contact_radius * a2 / gm->Reff -
      MY_2PI * a2 * sqrt(4.0 * cohesion * Emix / (MY_PI * gm->contact_radius));
  F_pulloff = 3.0 * MY_PI * cohesion * gm->Reff;

  return Fne;
}

/* ---------------------------------------------------------------------- */

void GranSubModNormalJKR::set_fncrit()
{
  Fncrit = fabs(Fne + 2.0 * F_pulloff);
}

/* ----------------------------------------------------------------------
   MDR contact model

   Contributing authors:
   William Zunker (MIT), Sachith Dunatunga (MIT),
   Dan Bolintineanu (SNL), Joel Clemmer (SNL)
------------------------------------------------------------------------- */

GranSubModNormalMDR::GranSubModNormalMDR(GranularModel *gm, LAMMPS *lmp) :
    GranSubModNormal(gm, lmp)
{
  num_coeffs = 6; // Young's Modulus, Poisson's ratio, yield stress, effective surface energy, psi_b, coefficent of restitution
  contact_radius_flag = 1;
  size_history = 39; 

  nondefault_history_transfer = 1;
  transfer_history_factor = new double[size_history];
  //transfer_history_factor[0] = +1;
  for (int i = 0; i < size_history; i++) { 
    transfer_history_factor[i] = +1;
  }
}

/* ---------------------------------------------------------------------- */

void GranSubModNormalMDR::coeffs_to_local()
{
  E = coeffs[0];      // Young's modulus
  nu = coeffs[1];     // Poisson's ratio
  Y = coeffs[2];      // yield stress
  gamma = coeffs[3];  // effective surface energy
  psi_b = coeffs[4];  // bulk response trigger based on ratio of remaining free area: A_{free}/A_{total}
  CoR = coeffs[5];    // coefficent of restitution

  if (E <= 0.0) error->all(FLERR, "Illegal MDR normal model, Young's modulus must be greater than 0");
  if (nu < 0.0 || nu > 0.5) error->all(FLERR, "Illegal MDR normal model, Poisson's ratio must be between 0 and 0.5");
  if (Y < 0.0) error->all(FLERR, "Illegal MDR normal model, yield stress must be greater than or equal to 0");
  if (gamma < 0.0) error->all(FLERR, "Illegal MDR normal model, effective surface energy must be greater than or equal to 0");
  if (psi_b < 0.0 || psi_b > 1.0) error->all(FLERR, "Illegal MDR normal model, psi_b must be between 0 and 1.0");
  if (CoR < 0.0 || CoR > 1.0) error->all(FLERR, "Illegal MDR normal model, coefficent of restitution must be between 0 and 1.0");
}

/* ---------------------------------------------------------------------- */

double GranSubModNormalMDR::calculate_forces()

{
  // To understand the structure of the overall code it is important to consider 
  // the following:
  //
  // The MDR contact model was developed by imagining individual particles being 
  // squished between a number of rigid flats (references below). To allow  
  // for many interacting particles, we extend the idea of isolated particles surrounded 
  // by rigid flats. In particular, we imagine placing rigid flats at the overlap 
  // midpoints between particles. The force is calculated seperately on both sides
  // of the contact assuming interaction with a rigid flat. The two forces are then 
  // averaged on either side of the contact to determine the final force. If the 
  // contact is between a particle and wall then only one force evaluation is required.
  //  
  // Zunker and Kamrin, 2024, Part I: https://doi.org/10.1016/j.jmps.2023.105492
  // Zunker and Kamrin, 2024, Part II: https://doi.org/10.1016/j.jmps.2023.105493

  const int i_true = gm->i;               // true i particle index
  const int j_true = gm->j;               // true j particle index
  const double radi_true = gm->radi;      // true i particle initial radius
  const double radj_true = gm->radj;      // true j particle initial radius
    
  double F;                               // average force 
  double F0;                              // force on contact side 0
  double F1;                              // force on contact side 1
  double R0;
  double R1;
  double delta = gm->delta;               // apparent overlap
  //if (gm->contact_type == PAIR) delta = gm->delta/2.0; // half displacement to imagine interaction with rigid flat 
  
  //std::cout << "Normal force is called for: " << i_true << ", " << j_true << std::endl;

  //std::cout << "Contact model has been entered " << gm->contact_type << ", " << PAIR << ", " << WALL << ", " << WALLREGION << ", " << gm->itype << ", " << gm->jtype << ", " << gm->delta << std::endl; 

  // initialize indexing in history array of different constact history variables 
  const int delta_offset_0 = 0;           // apparent overlap 
  const int delta_offset_1 = 1;           
  const int deltao_offset_0 = 2;          // displacement 
  const int deltao_offset_1 = 3;  
  const int delta_MDR_offset_0 = 4;       // MDR apparent overlap
  const int delta_MDR_offset_1 = 5;
  const int delta_BULK_offset_0 = 6;      // bulk displacement
  const int delta_BULK_offset_1 = 7;
  const int deltamax_MDR_offset_0 = 8;    // maximum MDR apparent overlap 
  const int deltamax_MDR_offset_1 = 9;
  const int Yflag_offset_0 = 10;          // yield flag
  const int Yflag_offset_1 = 11;
  const int deltaY_offset_0 = 12;         // yield displacement
  const int deltaY_offset_1 = 13;
  const int cA_offset_0 = 14;             // contact area intercept
  const int cA_offset_1 = 15;
  const int aAdh_offset_0 = 16;           // adhesive contact radius
  const int aAdh_offset_1 = 17;
  const int Ac_offset_0 = 18;             // contact area
  const int Ac_offset_1 = 19;
  const int eps_bar_offset_0 = 20;        // volume-averaged infinitesimal strain tensor
  const int eps_bar_offset_1 = 21;
  const int penalty_offset_ = 22;         // contact penalty 

  // temporary contact history variables
  const int deltamax_offset_ = 23;
  const int delta2_offset_ = 24;
  const int F_offset_0 = 25;
  const int F_offset_1 = 26;
  const int delta2_offset_0 = 27;           // apparent overlap 
  const int delta2_offset_1 = 28; 
  const int h_offset_0 = 29;
  const int h_offset_1 = 30;  
  const int deltae1D_offset_0 = 31;
  const int deltae1D_offset_1 = 32;  
  const int h_BULK_offset_0 = 33;
  const int h_BULK_offset_1 = 34;  
  const int k_BULK_offset_0 = 35;
  const int k_BULK_offset_1 = 36;
  const int k_MDR_offset_0 = 37;
  const int k_MDR_offset_1 = 38; 

  // initialize particle history variables 
  int tmp1, tmp2;
  int index_Ro = atom->find_custom("Ro",tmp1,tmp2);                       // initial radius
  int index_Vcaps = atom->find_custom("Vcaps",tmp1,tmp2);                 // spherical cap volume from intersection of apparent radius particle and contact planes
  int index_Vgeo = atom->find_custom("Vgeo",tmp1,tmp2);                   // geometric particle volume of apparent particle after removing spherical cap volume
  int index_Velas = atom->find_custom("Velas",tmp1,tmp2);                 // particle volume from linear elasticity  
  int index_eps_bar = atom->find_custom("eps_bar",tmp1,tmp2);             // volume-averaged infinitesimal strain tensor
  int index_dRnumerator = atom->find_custom("dRnumerator",tmp1,tmp2);     // summation of numerator terms in calculation of dR
  int index_dRdenominator = atom->find_custom("dRdenominator",tmp1,tmp2); // summation of denominator terms in calculation of dR
  int index_Acon0 = atom->find_custom("Acon0",tmp1,tmp2);                 // total area involved in contacts: Acon^{n} 
  int index_Acon1 = atom->find_custom("Acon1",tmp1,tmp2);                 // total area involved in contacts: Acon^{n+1}
  int index_Atot = atom->find_custom("Atot",tmp1,tmp2);                   // total particle area 
  int index_Atot_sum = atom->find_custom("Atot_sum",tmp1,tmp2);           // running sum of contact area minus cap area
  int index_ddelta_bar = atom->find_custom("ddelta_bar",tmp1,tmp2);       // change in mean surface displacement
  int index_psi = atom->find_custom("psi",tmp1,tmp2);                     // ratio of free surface area to total surface area
  int index_sigmaxx = atom->find_custom("sigmaxx",tmp1,tmp2);             // xx-component of the stress tensor, not necessary for force calculation
  int index_sigmayy = atom->find_custom("sigmayy",tmp1,tmp2);             // yy-component of the stress tensor, not necessary for force calculation  
  int index_sigmazz = atom->find_custom("sigmazz",tmp1,tmp2);             // zz-component of the stress tensor, not necessary for force calculation 
  int index_contacts = atom->find_custom("contacts",tmp1,tmp2);                     // total contacts on particle 
  int index_adhesive_length = atom->find_custom("adhesive_length",tmp1,tmp2);       // total contacts on particle 
  double * Rinitial = atom->dvector[index_Ro];
  double * Vgeo = atom->dvector[index_Vgeo];
  double * Velas = atom->dvector[index_Velas];
  double * Vcaps = atom->dvector[index_Vcaps];
  double * eps_bar = atom->dvector[index_eps_bar];
  double * dRnumerator = atom->dvector[index_dRnumerator];
  double * dRdenominator = atom->dvector[index_dRdenominator];
  double * Acon0 = atom->dvector[index_Acon0];
  double * Acon1 = atom->dvector[index_Acon1]; 
  double * Atot = atom->dvector[index_Atot];
  double * Atot_sum = atom->dvector[index_Atot_sum];
  double * ddelta_bar = atom->dvector[index_ddelta_bar];
  double * psi = atom->dvector[index_psi];
  double * sigmaxx = atom->dvector[index_sigmaxx];
  double * sigmayy = atom->dvector[index_sigmayy];
  double * sigmazz = atom->dvector[index_sigmazz];
  double * contacts = atom->dvector[index_contacts];
  double * adhesive_length = atom->dvector[index_adhesive_length];
  

  double * history = & gm->history[history_index]; // load in all history variables  

  // Rigid flat placement schemes
  double * Yflag_offset0 = & history[Yflag_offset_0];
  double * Yflag_offset1 = & history[Yflag_offset_1];
  double Yflag0 = *Yflag_offset0;
  double Yflag1 = *Yflag_offset1;
  double * cA_offset0 = & history[cA_offset_0];
  double * cA_offset1 = & history[cA_offset_1];
  double cA0 = *cA_offset0;
  double cA1 = *cA_offset1;

  // More rigid flat placement schemes
  double * a_offset0 = & history[aAdh_offset_0];
  double * a_offset1 = & history[aAdh_offset_1];
  double a0 = *a_offset0;
  double a1 = *a_offset1;
  double * F_offset0 = & history[F_offset_0];
  double * F_offset1 = & history[F_offset_1];
  double F0old = *F_offset0;
  double F1old = *F_offset1;
  double * h_offset0 = & history[h_offset_0];
  double * h_offset1 = & history[h_offset_1];
  double h0 = *h_offset0;
  double h1 = *h_offset1;
  double * h_BULK_offset0 = & history[h_BULK_offset_0];
  double * h_BULK_offset1 = & history[h_BULK_offset_1];
  double h_BULK0 = *h_BULK_offset0;
  double h_BULK1 = *h_BULK_offset1;
  double * k_BULK_offset0 = & history[k_BULK_offset_0];
  double * k_BULK_offset1 = & history[k_BULK_offset_1];
  double k_BULK0 = *k_BULK_offset0;
  double k_BULK1 = *k_BULK_offset1;
  double * delta2_offset0 = & history[delta2_offset_0];
  double * delta2_offset1 = & history[delta2_offset_1];
  double * delta2_offset = & history[delta2_offset_];
  double ddelta = gm->delta - *delta2_offset;
  *delta2_offset = gm->delta;
  double * k_MDR_offset0 = & history[k_MDR_offset_0];
  double * k_MDR_offset1 = & history[k_MDR_offset_1];
  double k0 = *k_MDR_offset0;
  double k1 = *k_MDR_offset1;
  *k_MDR_offset0 = 2*E/(1.0-pow(nu,2.0));
  *k_MDR_offset1 = 2*E/(1.0-pow(nu,2.0));
  double dde0;
  double dde1;
  if ((a0*h0 == 0.0 && k_BULK0*h_BULK0 == 0.0 && a1*h1 == 0.0 && k_BULK1*h_BULK1 == 0.0) || (F0old == 0.0 && F1old == 0.0)) {
    //std::cout << "Rigid flat placement case 1" << std::endl;
    dde0 = ddelta/2.0;
    dde1 = ddelta/2.0;
  } else if ((a0*h0 == 0.0 && k_BULK0*h_BULK0 == 0.0) || F0old == 0.0) {
    //std::cout << "Rigid flat placement case 2" << std::endl;
    dde0 = ddelta/2.0;
    dde1 = ddelta/2.0;
  } else if ((a1*h1 == 0.0 && k_BULK1*h_BULK1 == 0.0) || F1old == 0.0){
    //std::cout << "Rigid flat placement case 3" << std::endl;
    dde0 = ddelta/2.0;
    dde1 = ddelta/2.0;
  } else {
    //std::cout << "Rigid flat placement case 4" << std::endl;
    dde0 = -((F0old - F1old - a1*ddelta*k1*h1 - ddelta*h_BULK1*k_BULK1)/(a0*k0*h0 + h_BULK0*k_BULK0 + a1*k1*h1 + h_BULK1*k_BULK1));
    dde1 = -((F1old - F0old - a0*ddelta*k0*h0 - ddelta*h_BULK0*k_BULK0)/(a0*k0*h0 + h_BULK0*k_BULK0 + a1*k1*h1 + h_BULK1*k_BULK1));
  }
  //if (abs(dde0) > abs(ddelta) ||  abs(dde1) > abs(ddelta)) {
  //  dde0 = ddelta/2.0;
  //  dde1 = ddelta/2.0;
  //}

  //(*delta2_offset0 == 0.0 || a0 == 0.0 || a1 == 0.0) ? dde0 = ddelta/2.0 : dde0 = -((F0old - F1old - a1*ddelta*k1*h1 - ddelta*h_BULK1*k_BULK1)/(a0*k0*h0 + h_BULK0*k_BULK0 + a1*k1*h1 + h_BULK1*k_BULK1));
  //(*delta2_offset1 == 0.0 || a0 == 0.0 || a1 == 0.0) ? dde1 = ddelta/2.0 : dde1 = -((F1old - F0old - a0*ddelta*k0*h0 - ddelta*h_BULK0*k_BULK0)/(a0*k0*h0 + h_BULK0*k_BULK0 + a1*k1*h1 + h_BULK1*k_BULK1));
  //(*delta2_offset0 == 0.0) ? dde0 = ddelta/2.0 : dde0 = -((F0old - F1old - a1*ddelta*k1*h1)/(a0*k0*h0 + a1*k1*h1));
  //(*delta2_offset1 == 0.0) ? dde1 = ddelta/2.0 : dde1 = -((F1old - F0old  - a0*ddelta*k0*h0)/(a0*k0*h0 + a1*k1*h1));
  double delta0 = *delta2_offset0 + dde0;
  double delta1 = *delta2_offset1 + dde1;
  //if (i_true == 167 && j_true == 204) {
  //std::cout << "Normal model: " << gm->delta << ", " << ddelta << ", " << gm->radi << ", " << gm->radj << " | delta: " << delta0 << ", " << delta1 << " | delta2_offset: " << *delta2_offset0 << ", " << *delta2_offset1 << "| dde: " << dde0 << ", " << dde1 << ", " << (dde0 + dde1) << "| Fold: " << F0old << ", " << F1old << " | a: " << a0 << ", " << a1 << " | h: " << h0 << ", " << h1 << " | k_BULK: " << k_BULK0 << ", " << k_BULK1 << " | h_BULK: " << h_BULK0 << ", " << h_BULK1 << std::endl;
  //}
  *delta2_offset0 = delta0;
  *delta2_offset1 = delta1;

  

  double * deltamax_offset = & history[deltamax_offset_];
  if (gm->delta > *deltamax_offset) *deltamax_offset = gm->delta;
  double deltamax = *deltamax_offset;

  //std::cout << "Yflag: " << Yflag0 << ", " << Yflag1 << std::endl;

  for (int contactSide = 0; contactSide < 2; contactSide++) { 

    double * delta_offset; 
    double * deltao_offset;
    double * delta_MDR_offset;   
    double * delta_BULK_offset; 
    double * deltamax_MDR_offset; 
    double * Yflag_offset; 
    double * deltaY_offset; 
    double * cA_offset;
    double * aAdh_offset; 
    double * Ac_offset; 
    double * eps_bar_offset; 
    double * penalty_offset;

    // added for rigid flat placement
    double * h_offset;
    double * deltae1D_offset;
    double * h_BULK_offset;
    double * k_BULK_offset;

    if (contactSide == 0) {
      if (gm->contact_type == PAIR) {
        gm->i = std::max(i_true,j_true);
        gm->j = std::min(i_true,j_true);
        if (gm->i == i_true) {
          gm->radi = radi_true;
          gm->radj = radj_true;
        } else {
          gm->radi = radj_true;
          gm->radj = radi_true;
        }
        R0 = gm->radi;
        R1 = gm->radj;

        
        delta = delta0;

        //if (Yflag0 == 0.0 && Yflag1 == 0.0) {
        //  double deltaOpt1 = gm->delta*R0/(R0+R1);
        //  double deltaOpt2 = gm->delta*R1/(R0+R1);
        //  (gm->radi < gm->radj) ? delta = MAX(deltaOpt1,deltaOpt2) : delta = MIN(deltaOpt1,deltaOpt2); 
        //  std::cout << "CS 0, Case 1: " << R0 << ", " << R1 << "| " << gm->delta*R1/(R0+R1) << ", " << gm->delta*R0/(R0+R1) << std::endl;
        //} else if(Yflag0 == 1.0 && Yflag1 == 1.0) {
        //  double deltaOpt1 = (cA0 - cA1 + pow(gm->delta,2.0)*M_PI - 2.0*gm->delta*M_PI*R1)/(2.0*M_PI*(gm->delta - R0 - R1));
        //  double deltaOpt2 = (cA1 - cA0 + pow(gm->delta,2.0)*M_PI - 2.0*gm->delta*M_PI*R0)/(2.0*M_PI*(gm->delta - R0 - R1));
        //  (gm->radi < gm->radj) ? delta = MAX(deltaOpt1,deltaOpt2) : delta = MIN(deltaOpt1,deltaOpt2);
        //  std::cout << "CS 0, Case 2: " << R0 << ", " << R1 << "| " << cA0 << ", " << cA1 << " | " << deltaOpt1 << ", " << deltaOpt2 << std::endl;
        //} else if(Yflag0 == 0.0 && Yflag1 == 1.0) {
        //  double deltaOpt1 = 1.0/2.0*(R0 - sqrt(-4.0*(-(cA1/M_PI) + gm->delta*R0) + pow(-R0 - 2.0*R1,2.0)) + 2.0*R1);
        //  double deltaOpt2 = gm->delta - 1.0/2.0*(R0 - sqrt(-4.0*(-(cA1/M_PI) + gm->delta*R0) + pow(-R0 - 2.0*R1,2.0)) + 2.0*R1);
        //  (gm->radi < gm->radj) ? delta = MAX(deltaOpt1,deltaOpt2) : delta = MIN(deltaOpt1,deltaOpt2);
        //  std::cout << "CS 0, Case 3: " << gm->radi << ", " << gm->radj << "| " << deltaOpt1 << ", " << deltaOpt2 << std::endl;
        //} else {
        //  double deltaOpt1 = 1.0/2.0*(2.0*gm->delta - 2.0*R0 - R1 + sqrt(4.0*cA0 + 4.0*M_PI*pow(R0,2.0) - 4.0*gm->delta*M_PI*R1 + 4.0*M_PI*R0*R1 + M_PI*pow(R1,2.0))/sqrt(M_PI));
        //  double deltaOpt2 = gm->delta - 1.0/2.0*(2.0*gm->delta - 2.0*R0 - R1 + sqrt(4.0*cA0 + 4.0*M_PI*pow(R0,2.0) - 4.0*gm->delta*M_PI*R1 + 4.0*M_PI*R0*R1 + M_PI*pow(R1,2.0))/sqrt(M_PI));
        //  (gm->radi < gm->radj) ? delta = MAX(deltaOpt1,deltaOpt2) : delta = MIN(deltaOpt1,deltaOpt2);
        //  std::cout << "CS 0, Case 4: " << gm->radi << ", " << gm->radj << "| " << deltaOpt1 << ", " << deltaOpt2 << std::endl;
        //}
        //std::cout << "Contact side 0: " << gm->radi << ", " << gm->radj << "| " << gm->delta << ", " << delta << std::endl;

        //double deltaOpt1 = gm->delta*(gm->delta - 2.0*R1)/(2.0*(gm->delta - R0 - R1));
        //double deltaOpt2 = gm->delta*(gm->delta - 2.0*R0)/(2.0*(gm->delta - R0 - R1));
        //(gm->radi < gm->radj) ? delta = MAX(deltaOpt1,deltaOpt2) : delta = MIN(deltaOpt1,deltaOpt2);

        //double deltaOpt1 = deltamax*(deltamax - 2.0*R1)/(2.0*(deltamax - R0 - R1));
        //double deltaOpt2 = deltamax*(deltamax - 2.0*R0)/(2.0*(deltamax - R0 - R1));
        //double delta_ratio;
        //(gm->radi < gm->radj) ? delta_ratio = MAX(deltaOpt1,deltaOpt2)/deltamax : delta_ratio = MIN(deltaOpt1,deltaOpt2)/deltamax;
        //delta = gm->delta*delta_ratio;

      }
      delta_offset = & history[delta_offset_0];
      deltao_offset = & history[deltao_offset_0];
      delta_MDR_offset = & history[delta_MDR_offset_0];
      delta_BULK_offset = & history[delta_BULK_offset_0];
      deltamax_MDR_offset = & history[deltamax_MDR_offset_0];
      Yflag_offset = & history[Yflag_offset_0];
      deltaY_offset = & history[deltaY_offset_0];
      cA_offset = & history[cA_offset_0];
      aAdh_offset = & history[aAdh_offset_0];
      Ac_offset = & history[Ac_offset_0];
      eps_bar_offset = & history[eps_bar_offset_0];

      // added for rigid flat placement
      h_offset = & history[h_offset_0];
      deltae1D_offset = & history[deltae1D_offset_0];
      h_BULK_offset = & history[h_BULK_offset_0];
      k_BULK_offset = & history[k_BULK_offset_0];
    } else {
      if (gm->contact_type != PAIR) break; // contact with particle-wall requires only one evaluation
      gm->i = std::min(i_true,j_true);
      gm->j = std::max(i_true,j_true);
      if (gm->i == i_true) {
        gm->radi = radi_true;
        gm->radj = radj_true;
      } else {
        gm->radi = radj_true;
        gm->radj = radi_true;
      }

      delta = delta1;

      //if (Yflag0 == 0.0 && Yflag1 == 0.0) {
      //  double deltaOpt1 = gm->delta*R0/(R0+R1);
      //  double deltaOpt2 = gm->delta*R1/(R0+R1);
      //  (gm->radi < gm->radj) ? delta = MAX(deltaOpt1,deltaOpt2) : delta = MIN(deltaOpt1,deltaOpt2);  
      //  std::cout << "CS 1, Case 1: " << R0 << ", " << R1 << "| " << gm->delta*R1/(R0+R1) << ", " << gm->delta*R0/(R0+R1) << std::endl;
      //} else if(Yflag0 == 1.0 && Yflag1 == 1.0) {
      //  double deltaOpt1 = (cA0 - cA1 + pow(gm->delta,2.0)*M_PI - 2.0*gm->delta*M_PI*R1)/(2.0*M_PI*(gm->delta - R0 - R1));
      //  double deltaOpt2 = (cA1 - cA0 + pow(gm->delta,2.0)*M_PI - 2.0*gm->delta*M_PI*R0)/(2.0*M_PI*(gm->delta - R0 - R1));
      //  (gm->radi < gm->radj) ? delta = MAX(deltaOpt1,deltaOpt2) : delta = MIN(deltaOpt1,deltaOpt2);
      //  std::cout << "CS 1, Case 2: " << R0 << ", " << R1 << "| " << cA0 << ", " << cA1 << " | " << deltaOpt1 << ", " << deltaOpt2 << std::endl;
      //} else if(Yflag0 == 0.0 && Yflag1 == 1.0) {
      //  double deltaOpt1 = 1.0/2.0*(R0 - sqrt(-4.0*(-(cA1/M_PI) + gm->delta*R0) + pow(-R0 - 2.0*R1,2.0)) + 2.0*R1);
      //  double deltaOpt2 = gm->delta - 1.0/2.0*(R0 - sqrt(-4.0*(-(cA1/M_PI) + gm->delta*R0) + pow(-R0 - 2.0*R1,2.0)) + 2.0*R1);
      //  (gm->radi < gm->radj) ? delta = MAX(deltaOpt1,deltaOpt2) : delta = MIN(deltaOpt1,deltaOpt2);
      //  std::cout << "CS 1, Case 3: " << gm->radi << ", " << gm->radj << "| " << deltaOpt1 << ", " << deltaOpt2 << std::endl;
      //} else {
      //  double deltaOpt1 = 1.0/2.0*(2.0*gm->delta - 2.0*R0 - R1 + sqrt(4.0*cA0 + 4.0*M_PI*pow(R0,2.0) - 4.0*gm->delta*M_PI*R1 + 4.0*M_PI*R0*R1 + M_PI*pow(R1,2.0))/sqrt(M_PI));
      //  double deltaOpt2 = gm->delta - 1.0/2.0*(2.0*gm->delta - 2.0*R0 - R1 + sqrt(4.0*cA0 + 4.0*M_PI*pow(R0,2.0) - 4.0*gm->delta*M_PI*R1 + 4.0*M_PI*R0*R1 + M_PI*pow(R1,2.0))/sqrt(M_PI));
      //  (gm->radi < gm->radj) ? delta = MAX(deltaOpt1,deltaOpt2) : delta = MIN(deltaOpt1,deltaOpt2);
      //  std::cout << "CS 1, Case 4: " << gm->radi << ", " << gm->radj << "| " << deltaOpt1 << ", " << deltaOpt2 << std::endl;
      //}

      //double deltaOpt1 = gm->delta*(gm->delta - 2.0*R1)/(2.0*(gm->delta - R0 - R1));
      //double deltaOpt2 = gm->delta*(gm->delta - 2.0*R0)/(2.0*(gm->delta - R0 - R1));
      //(gm->radi < gm->radj) ? delta = MAX(deltaOpt1,deltaOpt2) : delta = MIN(deltaOpt1,deltaOpt2);

      //double deltaOpt1 = deltamax*(deltamax - 2.0*R1)/(2.0*(deltamax - R0 - R1));
      //double deltaOpt2 = deltamax*(deltamax - 2.0*R0)/(2.0*(deltamax - R0 - R1));
      //double delta_ratio;
      //(gm->radi < gm->radj) ? delta_ratio = MAX(deltaOpt1,deltaOpt2)/deltamax : delta_ratio = MIN(deltaOpt1,deltaOpt2)/deltamax;
      //delta = gm->delta*delta_ratio;
      
      //delta = gm->delta*(gm->delta - 2.0*gm->radj)/(2.0*(gm->delta - gm->radj - gm->radi));
      //std::cout << "Contact side 1: " << gm->radi << ", " << gm->radj << "| " << gm->delta << ", " << delta << std::endl;
      delta_offset = & history[delta_offset_1];
      deltao_offset = & history[deltao_offset_1];
      delta_MDR_offset = & history[delta_MDR_offset_1];
      delta_BULK_offset = & history[delta_BULK_offset_1];
      deltamax_MDR_offset = & history[deltamax_MDR_offset_1];
      Yflag_offset = & history[Yflag_offset_1];
      deltaY_offset = & history[deltaY_offset_1];
      cA_offset = & history[cA_offset_1];
      aAdh_offset = & history[aAdh_offset_1];
      Ac_offset = & history[Ac_offset_1];
      eps_bar_offset = & history[eps_bar_offset_1];

      // added for rigid flat placement
      h_offset = & history[h_offset_1];
      deltae1D_offset = & history[deltae1D_offset_1];
      h_BULK_offset = & history[h_BULK_offset_1];
      k_BULK_offset = & history[k_BULK_offset_1];
    }

    //delta = gm->delta;

    // temporary i and j indices
    const int i = gm->i;
    const int j = gm->j;


    //if (lmp->update->ntimestep == 229999 && i == 2 && j == 0) {
    //  std::cout << "timestep is equal to 229999" << std::endl;
    //}

    //if ( (i == 0 && j == 2 && gm->contact_type == 0) || (i == 2 && j == 0 && gm->contact_type == 0)) {
    //if ( (gm->contact_type == 0) ) {
      //std::cout << i << ", " << j << ", " << gm->contact_type << " || " << *delta_offset << ", " << (uintptr_t)(delta_offset) << " || " << *deltao_offset << ", " << (uintptr_t)(deltao_offset) << " || " << *delta_MDR_offset << ", " << (uintptr_t)(delta_MDR_offset) << " || " << *delta_BULK_offset << ", " << (uintptr_t)(delta_BULK_offset) << " || " << *deltamax_MDR_offset << ", " << (uintptr_t)(deltamax_MDR_offset) << " || "  << *Yflag_offset << ", " << (uintptr_t)(Yflag_offset) << " || " << *deltaY_offset << ", " << (uintptr_t)(deltaY_offset) << " || " << *cA_offset << ", " << (uintptr_t)(cA_offset) << " || " << *aAdh_offset << ", " << (uintptr_t)(aAdh_offset) << " || " << *Ac_offset << ", " << (uintptr_t)(Ac_offset) << " || " << *eps_bar_offset << ", " << (uintptr_t)(eps_bar_offset) << std::endl;
    //  std::cout << i << ", " << j << ", " << gm->contact_type << " || " << *delta_offset << ", " << (uintptr_t)(delta_offset) << " || " << *deltao_offset << ", " << (uintptr_t)(deltao_offset) << " || " << *delta_MDR_offset << ", " << (uintptr_t)(delta_MDR_offset) << " || " << *delta_BULK_offset << ", " << (uintptr_t)(delta_BULK_offset) << " || " << *deltamax_MDR_offset << ", " << (uintptr_t)(deltamax_MDR_offset) << " || "  << *Yflag_offset << ", " << (uintptr_t)(Yflag_offset) << " || " << *deltaY_offset << ", " << (uintptr_t)(deltaY_offset) << " || " << Velas[i] << ", " << eps_bar[i] << std::endl;
    //}

    // material and geometric property definitions
    // E, nu, Y gamma , psi_b, and CoR are already defined.
    const double G = E/(2.0*(1.0+nu));          // shear modulus
    const double kappa = E/(3.0*(1.0-2.0*nu));  // bulk modulus
    const double Eeff = E/(1.0-pow(nu,2.0));    // composite plane strain modulus

    const double Ro = Rinitial[i];              // initial radius
    const double R = gm->radi;                  // apparent radius
 
    // kinematics 
    const double ddelta = delta - *delta_offset;
    *delta_offset = delta;

    const double deltao = delta - (R - Ro);
    const double ddeltao = deltao - *deltao_offset;
    *deltao_offset = deltao;

    double ddelta_MDR;
    double ddelta_BULK;
    if ( psi[i] < psi_b ) { // if true, bulk response has triggered, split displacement increment between the MDR and BULK components 
      ddelta_MDR = std::min(ddelta-ddelta_bar[i], delta-*delta_MDR_offset);
      ddelta_BULK = ddelta_bar[i];
    } else { // if false, no bulk response, full displacement increment goes to the MDR component
      ddelta_BULK = 0.0;
      ddelta_MDR = ddelta;
    }
    const double delta_MDR = *delta_MDR_offset + ddelta_MDR; // MDR displacement
    *delta_MDR_offset = delta_MDR; // Update old MDR displacement
    const double delta_BULK = std::max(0.0,*delta_BULK_offset+ddelta_BULK); // bulk displacement
    *delta_BULK_offset = delta_BULK; // update old bulk displacement

    if (delta_MDR > *deltamax_MDR_offset) *deltamax_MDR_offset = delta_MDR;
    const double deltamax_MDR = *deltamax_MDR_offset;

    const double pY = Y*(1.75*exp(-4.4*deltamax_MDR/R) + 1.0); // Set value of average pressure along yield surface
    if ( *Yflag_offset == 0.0 && delta_MDR >= deltamax_MDR ) {
    const double phertz = 4*Eeff*sqrt(delta_MDR)/(3*M_PI*sqrt(R));
      if ( phertz > pY ) {
        *Yflag_offset = 1.0;
        *deltaY_offset = delta_MDR;
        *cA_offset = M_PI*(pow(*deltaY_offset,2.0) - *deltaY_offset*R);
      }
    } 

    //if (i_true == 167 && j_true == 204) {
    //std::cout << "i " << i << " | j " << j << " | delta_BULK: " << delta_BULK << " | delta_MDR " << delta_MDR << " | ddelta_BULK " << ddelta_BULK << " | ddelta_MDR " << ddelta_MDR << std::endl;
    //}

    //std::cout << "Yield Flag: " << *Yflag_offset << ", " << R << std::endl;

    // MDR force calculation
    double F_MDR;
    double A;                     // height of elliptical indenter
    double B;                     // width of elliptical indenter
    double deltae1D;              // transformed elastic displacement
    double deltaR;                // displacement correction 
    double amax;                  // maximum experienced contact radius
    const double cA = *cA_offset; // contact area intercept

    if ( *Yflag_offset == 0.0 ) { // elastic contact
      A = 4.0*R;              
      B = 2.0*R;              
      deltae1D = delta_MDR;    
      (deltae1D > 0) ?  amax = sqrt(deltae1D*R) : amax = 0.0;  
    } else { // plastic contact
      amax = sqrt((2.0*deltamax_MDR*R - pow(deltamax_MDR,2.0)) + cA/M_PI);              
      A = 4.0*pY/Eeff*amax;                                                             
      B = 2.0*amax;                                                                     
      const double deltae1Dmax = A/2.0; // maximum transformed elastic displacement 
      const double Fmax = Eeff*(A*B/4.0)*(acos(1.0 - 2.0*deltae1Dmax/A) - (1.0 - 2.0*deltae1Dmax/A)*sqrt(4.0*deltae1Dmax/A - 4.0*pow(deltae1Dmax,2.0)/pow(A,2.0))); // force caused by full submersion of elliptical indenter to depth of A/2
      const double zR = R - (deltamax_MDR - deltae1Dmax); // depth of particle center
      deltaR = (Fmax*(2*pow(amax,2.0)*(-1 + nu) - (-1 + 2*nu)*zR*(-zR + sqrt(pow(amax,2.0) + pow(zR,2.0)))))/((M_PI*pow(amax,2.0))*2*G*sqrt(pow(amax,2.0) + pow(zR,2.0))); 
      deltae1D = (delta_MDR - deltamax_MDR + deltae1Dmax + deltaR)/(1 + deltaR/deltae1Dmax);  // transformed elastic displacement 
    }

    double ddeltae1D = deltae1D - *deltae1D_offset;
    *deltae1D_offset = deltae1D;
    (ddelta != 0.0 ) ? *h_offset = abs(ddeltae1D/ddelta) : *h_offset = 1.0;
    //if (i_true == 4 && j_true == 52){
    //std::cout << "h_offset: " << *h_offset << " | ddeltae1D: " << ddeltae1D << " | ddelta: " << ddelta << " | R: " << R << std::endl;
    //}

    //std::cout << psi_b << ", " << psi[i] << ", " << A << ", " << B << ", " << pY << ", " << amax << " || " << deltao << ", " << delta << ", " << ddelta << ", " << *delta_offset << ", " << ddelta_bar[i] << " || " << delta_MDR << ", " << ddelta_MDR << ", " << *delta_MDR_offset << ", " << deltamax_MDR << " || " << delta_BULK << ", " << ddelta_BULK << ", " << *delta_BULK_offset << " || " << R << std::endl;
    
    //std::cout << i << ", " << j << ", " << A << ", " << B << " || " << deltao << ", " << delta << ", " << ddelta << ", " << R <<  ", " << M_PI*pow(amax,2.0) << std::endl;

    double a_na;
    (deltae1D >= 0.0) ? a_na = B*sqrt(A - deltae1D)*sqrt(deltae1D)/A : a_na = 0.0;
    double aAdh = *aAdh_offset; 

    //if (i_true == 4 && j_true == 52){
    //std::cout << "CS: " << contactSide << ", aAdh: " << aAdh << ", deltae1D: " << deltae1D << ", A: " << A << ", B:" << B << ", amax: " << amax << ", deltae1D: " << deltae1D << ", R: " << R << std::endl;
    //}


    if ( gamma > 0.0  ) { // adhesive contact
    double g_aAdh;
    
      if (delta_MDR == deltamax_MDR || a_na >= aAdh ) { // case 1: no tensile springs, purely compressive contact
        (deltae1D <= 0.0) ? F_MDR = 0.0 : F_MDR = Eeff*(A*B/4.0)*(acos(1.0 - 2.0*deltae1D/A) - (1.0 - 2.0*deltae1D/A)*sqrt(4.0*deltae1D/A - 4.0*pow(deltae1D,2.0)/pow(A,2.0))); 
        if ( std::isnan(F_MDR) ) {
           std::cout << "F_MDR is NaN, case 1: no tensile springs" << std::endl;
           std::cout << "Normal model: " << gm->delta << ", " << ddelta << ", " << gm->radi << ", " << gm->radj << " | delta: " << delta0 << ", " << delta1 << " | delta2_offset: " << *delta2_offset0 << ", " << *delta2_offset1 << "| dde: " << dde0 << ", " << dde1 << "| Fold: " << F0old << ", " << F1old << " | a: " << a0 << ", " << a1 << " | k_BULK: " << k_BULK0 << ", " << k_BULK1 << " | h_BULK: " << h_BULK0 << ", " << h_BULK1 << std::endl;
           std::cout << "i_true: " << i_true << "j_true: " << j_true << "deltae1D: " << deltae1D << ", A: " << A << ", B:" << B << ", amax:" << amax << ", deltamax_MDR" << deltamax_MDR << ", R:" << R << std::endl;
           std::exit(1);
        }
        *aAdh_offset = a_na;
      } else {
        const double lmax = sqrt(2.0*M_PI*aAdh*gamma/Eeff); 
        g_aAdh = A/2 - A/B*sqrt(pow(B,2.0)/4 - pow(aAdh,2.0));
        const double acrit = (-((pow(B,2)*gamma*M_PI)/(pow(A,2)*Eeff)) + (pow(2,0.3333333333333333)*pow(B,4)*pow(gamma,2)*pow(M_PI,1.6666666666666667))/
                            (pow(A,2)*pow(Eeff,2)*pow((27*pow(A,4)*pow(B,4)*gamma)/Eeff - (2*pow(B,6)*pow(gamma,3)*pow(M_PI,2))/pow(Eeff,3) + (3*sqrt(3)*sqrt(27*pow(A,8)*pow(B,8)*pow(Eeff,2)*pow(gamma,2) -
                              4*pow(A,4)*pow(B,10)*pow(gamma,4)*pow(M_PI,2)))/pow(Eeff,2),0.3333333333333333)) + (pow(M_PI/2.,0.3333333333333333)*pow((27*pow(A,4)*pow(B,4)*gamma)/Eeff -
                            (2*pow(B,6)*pow(gamma,3)*pow(M_PI,2))/pow(Eeff,3) + (3*sqrt(3)*sqrt(27*pow(A,8)*pow(B,8)*pow(Eeff,2)*pow(gamma,2) - 4*pow(A,4)*pow(B,10)*pow(gamma,4)*pow(M_PI,2)))/
                            pow(Eeff,2),0.3333333333333333))/pow(A,2))/6;

        if ( (deltae1D + lmax - g_aAdh) >= 0.0) { // case 2: tensile springs, but not exceeding critical length --> deltae + lmax - g(aAdhes) >= 0
        //if (contactSide == 0) {
          //std::cout << "Case 2 tensile springs not exceeding critical length, R " << R <<  "deltae1D " << deltae1D << " , lmax " << lmax << " , g_adh" << g_aAdh << " sum " << (deltae1D + lmax - g_aAdh) << std::endl; 
        //}
          const double deltaeAdh = g_aAdh; 
          const double F_na = Eeff*(A*B/4.0)*(acos(1.0 - 2.0*deltaeAdh/A) - (1.0 - 2.0*deltaeAdh/A)*sqrt(4.0*deltaeAdh/A - 4.0*pow(deltaeAdh,2.0)/pow(A,2.0)));
          const double F_Adhes = 2.0*Eeff*(deltae1D - deltaeAdh)*aAdh;
          F_MDR = F_na + F_Adhes; 
          if ( std::isnan(F_MDR) ) std::cout << "F_MDR is NaN, case 2: tensile springs, but not exceeding critical length" << std::endl;

        } else { // case 3: tensile springs exceed critical length --> deltae + lmax - g(aAdhes) = 0
          //if (contactSide == 0) {
            //std::cout << "Case 3 tensile springs exceed critical length, R " << R << " deltae1D " << deltae1D << " , lmax " << lmax << " , g_adh " << g_aAdh << " sum " << (deltae1D + lmax - g_aAdh) << std::endl; 
          //}
        
          if ( aAdh < acrit ) {
            F_MDR = 0.0;
          } else {

            // bisection to find aAdh

            const double maxIterations = 100;
            const double error = 1e-8;
            double a_bisec = aAdh;
            double b_bisec = 0.9*acrit;
            double root = (a_bisec + b_bisec)/2.0;
            double froot;
            double fb_bisec;
            //if (aAdh = amax) {
            //  aAdh = aAdh - (aAdh-acrit)*0.001;
            //} else {
              for (int lv1 = 0; lv1 < maxIterations; ++lv1) {
                froot = deltae1D + sqrt(2.0*M_PI*root*gamma/Eeff) - ( A/2 - A/B*sqrt(pow(B,2.0)/4 - pow(root,2.0)) );
                if (abs(froot) < error) {
                  break;
                } else {
                  fb_bisec = deltae1D + sqrt(2.0*M_PI*b_bisec*gamma/Eeff) - ( A/2 - A/B*sqrt(pow(B,2.0)/4 - pow(b_bisec,2.0)) );
                  if (froot > 0.0 && fb_bisec > 0.0) {
                    b_bisec = root;
                  } else {
                    a_bisec = root;
                  }
                  root = (a_bisec + b_bisec)/2.0;
                }
                if (i == 123  && j == 166) {
                  std::cout << "Adhesion calculation - i " << i << ", j " << j << ", lv1 " << lv1 << ", froot " << froot << ", fb_bisec " << fb_bisec << ", root " << root << ", amax " << amax << ", acrit " << acrit << ", a_bisec " << a_bisec << ", b_bisec " << b_bisec << std::endl;
                }
                if (lv1 == maxIterations-1){
                  aAdh = 0.0;
                  std::cout << "Max iterations reached - i " << i << ", j " << j << ", froot " << froot << ", root " << root << ", amax " << amax << ", acrit: " << acrit << std::endl;
                  //std::exit(1);
                }
              }
            //}
            (root < acrit ) ? aAdh = 0.0 : aAdh = root;


            // newton-raphson to find aAdh

            //// newton-raphson to find aAdh
            //if (aAdh = amax) {
            //  aAdh = aAdh - (aAdh-acrit)*0.01;
            //} else {
            //  const double maxIterations = 100;
            //  const double error = 1e-10;
            //  const double error2 = 1e-16;
            //  double aAdh_tmp = aAdh;
            //  double fa; 
            //  double fa2;
            //  double dfda;
            //  for (int lv1 = 0; lv1 < maxIterations; ++lv1) {
            //    fa = deltae1D + sqrt(2.0*M_PI*aAdh_tmp*gamma/Eeff) - ( A/2 - A/B*sqrt(pow(B,2.0)/4 - pow(aAdh_tmp,2.0)) );
            //    if (abs(fa) < error) {
            //      break;
            //    } 
            //    dfda = -((aAdh_tmp*A)/(B*sqrt(-pow(aAdh_tmp,2.0) + pow(B,2.0)/4.0))) + (gamma*sqrt(M_PI/2.0))/(Eeff*sqrt((aAdh_tmp*gamma)/Eeff));
            //    aAdh_tmp = aAdh_tmp - fa/dfda;
            //    fa2 = deltae1D + sqrt(2.0*M_PI*aAdh_tmp*gamma/Eeff) - ( A/2 - A/B*sqrt(pow(B,2.0)/4 - pow(aAdh_tmp,2.0)) );
            //    if (abs(fa-fa2) < error2) {
            //      break;
            //    } 
            //    if (lv1 == maxIterations-1){
            //      aAdh_tmp = 0.0;
            //    }
            //  }
            //  aAdh = aAdh_tmp;
            //}

            //const double maxIterations = 1000;
            //const double error = 1e-8;
            //const double error2 = 1e-16;
            //double aAdh_tmp = aAdh;
            //double fa; 
            //double fa2;
            //double dfda;
            //if (aAdh = amax) {
            //  aAdh = aAdh - (aAdh-acrit)*0.01;
            //} else {
            //  for (int lv1 = 0; lv1 < maxIterations; ++lv1) {
            //    fa = deltae1D + sqrt(2.0*M_PI*aAdh_tmp*gamma/Eeff) - ( A/2 - A/B*sqrt(pow(B,2.0)/4 - pow(aAdh_tmp,2.0)) );
            //    if (abs(fa) < error) {
            //      //std::cout << "abs(fa) < error, fa " << fa << " CS " << contactSide << " lv1 " << lv1 << " R " << R << " deltae1D " << deltae1D << " , lmax " << sqrt(2.0*M_PI*aAdh_tmp*gamma/Eeff) << " , g_adh " << ( A/2 - A/B*sqrt(pow(B,2.0)/4 - pow(aAdh_tmp,2.0)) ) << std::endl;
            //      break;
            //    } 
            //    const double dfda_term1_denom = -pow(aAdh_tmp,2.0) + pow(B,2.0)/4.0;
            //    if (dfda_term1_denom == 0.0) {
            //      const double facrit = deltae1D + sqrt(2.0*M_PI*acrit*gamma/Eeff) - ( A/2 - A/B*sqrt(pow(B,2.0)/4 - pow(acrit,2.0)) ); 
            //      dfda = (facrit-fa)/(acrit-aAdh_tmp);
            //      std::cout << " CS " << contactSide << "facrit " << facrit << " fa " << fa << " acrit " << acrit  << " aAdh_tmp " << aAdh_tmp << " dfda " << dfda << std::endl;
            //    } else {
            //      dfda = -((aAdh_tmp*A)/(B*sqrt(-pow(aAdh_tmp,2.0) + pow(B,2.0)/4.0))) + (gamma*sqrt(M_PI/2.0))/(Eeff*sqrt((aAdh_tmp*gamma)/Eeff));
            //    }
            //    if ( (aAdh_tmp - fa/dfda < B/2.0) && (aAdh_tmp - fa/dfda > 0.0) ) aAdh_tmp = aAdh_tmp - fa/dfda;
            //    fa2 = deltae1D + sqrt(2.0*M_PI*aAdh_tmp*gamma/Eeff) - ( A/2 - A/B*sqrt(pow(B,2.0)/4 - pow(aAdh_tmp,2.0)) );
            //    if (abs(fa-fa2) < error2) {
            //      //std::cout << "abs(fa-fa2) < error, fa " << fa << " fa2 " << fa2 << " dfda " << dfda << " CS " << contactSide << " lv1 " << lv1 << " R " << R << " deltae1D " << deltae1D << " , lmax " << sqrt(2.0*M_PI*aAdh_tmp*gamma/Eeff) << " , g_adh " << ( A/2 - A/B*sqrt(pow(B,2.0)/4 - pow(aAdh_tmp,2.0)) ) << std::endl;
            //      //std::cout << "dfda terms: term 1 " << -((aAdh_tmp*A)/(B*sqrt(-pow(aAdh_tmp,2.0) + pow(B,2.0)/4.0))) << " term 2 " << (gamma*sqrt(M_PI/2.0))/(Eeff*sqrt((aAdh_tmp*gamma)/Eeff)) << " denom " << (-pow(aAdh_tmp,2.0) + pow(B,2.0)/4.0) << std::endl;
            //      break;
            //    } 
            //    if (lv1 == maxIterations-1){
            //      aAdh_tmp = 0.0;
            //      std::cout << "Max iterations reached: fa " << fa << " aAdh_tmp, " << aAdh_tmp << std::endl;
            //    }
            //  }
            //  aAdh = aAdh_tmp;
            //}

            //std::cout << " CS " << contactSide << " R " << R << "deltae1D " << deltae1D << " , lmax " << lmax << " , g_adh" << g_aAdh << ", aAdh " << aAdh << std::endl; 


          
            g_aAdh = A/2.0 - A/B*sqrt(pow(B,2.0)/4.0 - pow(aAdh,2.0));                  
            const double deltaeAdh = g_aAdh; 
            const double F_na = Eeff*(A*B/4.0)*(acos(1.0 - 2.0*deltaeAdh/A) - (1.0 - 2.0*deltaeAdh/A)*sqrt(4.0*deltaeAdh/A - 4.0*pow(deltaeAdh,2.0)/pow(A,2.0)));
            const double F_Adhes = 2.0*Eeff*(deltae1D - deltaeAdh)*aAdh;
            F_MDR = F_na + F_Adhes; 
            if ( std::isnan(F_MDR) ) std::cout << "F_MDR is NaN, case 3: tensile springs exceed critical length" << std::endl;
          }
          *aAdh_offset = aAdh;
        }
      }
    } else { // non-adhesive contact
      *aAdh_offset = a_na;
      (deltae1D <= 0.0) ? F_MDR = 0.0 : F_MDR = Eeff*(A*B/4.0)*(acos(1.0 - 2.0*deltae1D/A) - (1.0 - 2.0*deltae1D/A)*sqrt(4.0*deltae1D/A - 4.0*pow(deltae1D,2.0)/pow(A,2.0))); 
      if ( std::isnan(F_MDR) ) {
        std::cout << "F_MDR is NaN, non-adhesive case" << std::endl;
        std::cout << "Normal model: " << gm->delta << ", " << ddelta << ", " << gm->radi << ", " << gm->radj << " | delta: " << delta0 << ", " << delta1 << " | delta2_offset: " << *delta2_offset0 << ", " << *delta2_offset1 << "| dde: " << dde0 << ", " << dde1 << "| Fold: " << F0old << ", " << F1old << " | a: " << a0 << ", " << a1 << " | k_BULK: " << k_BULK0 << ", " << k_BULK1 << " | h_BULK: " << h_BULK0 << ", " << h_BULK1 << std::endl;
        std::cout << "i_true: " << i_true << "j_true: " << j_true << "deltae1D: " << deltae1D << ", A: " << A << ", B:" << B << ", amax:" << amax << ", deltamax_MDR" << deltamax_MDR << ", R:" << R << std::endl;
        std::exit(1);
      } 
    }
    
    //std::cout << gm->i << ", " << gm->j << ", aAdh_offset: " << *aAdh_offset << ", aAdh: " << aAdh << ", a_na: " << a_na << std::endl;

    // added for rigid flat placement
    (contactSide == 0) ? a0 = a_na : a1 = a_na;

    contacts[i] += 1;
    adhesive_length[i] += aAdh;

    // contact penalty scheme
    penalty_offset = & history[penalty_offset_];
    double pij = *penalty_offset;
    const double wij = std::max(1.0-pij,0.0);

    //std::cout << gm->i << ", " << gm->j << ", " << gm->contact_type << ", " << *gm->xi[1] << ", " << *gm->xj[1] << std::endl;

    // area related calculations 
    double Ac; 
    (*Yflag_offset == 0.0) ? Ac = M_PI*delta*R : Ac = M_PI*((2.0*delta*R - pow(delta,2.0)) + cA/M_PI);
    if (Ac < 0.0 ) Ac = 0.0;
    Atot_sum[i] += wij*(Ac - 2.0*M_PI*R*(deltamax_MDR + delta_BULK));
    Acon1[i] += wij*Ac;

    // bulk force calculation
    double F_BULK;
    (delta_BULK <= 0.0) ? F_BULK = 0.0 : F_BULK = (1.0/Vgeo[i])*Acon0[i]*delta_BULK*kappa*Ac;

    // added for rigid flat placement
    (Vgeo[i] == 0.0) ? *k_BULK_offset = 0.0 : *k_BULK_offset = (1.0/Vgeo[i])*Acon0[i]*kappa*Ac;
    (ddelta != 0.0) ? *h_BULK_offset = abs(ddelta_BULK/ddelta) : *h_BULK_offset = 1.0;
    //std::cout << "h_BULK_offset: " << *h_BULK_offset << " | ddeltae1D: " << ddelta_BULK << " | ddelta: " << ddelta << " | R: " << R << std::endl;
    //std::cout << "k_BULK_offset: " << *k_BULK_offset << ", Vgeo: " << Vgeo[i] << ", Acon0: " << Acon0[i] << ", kappa: " << kappa << ", Ac: " << Ac << ",wij :" << wij << std::endl;

    //if (F_BULK > 0.0) {
    //  std::cout << "F_BULK is: " << F_BULK << std::endl;
    //}

    //std::cout << delta_BULK << ", " << F_BULK << ", " << (1.0/Vgeo[i])*Acon0[i]*delta_BULK*kappa*Ac << std::endl;

    //std::cout << gm->i << ", " << gm->j << ", " << Vgeo[i] << ", " << Acon0[i] << ", " << Acon1[i] << ", " << Ac << ", " << kappa << " || " << psi[i] << ", " << ddelta_bar[i] << ", " << ddelta << ", " << ddelta_MDR << ", " << ddelta_BULK << ", " << delta << ", " << delta_MDR << ", " << delta_BULK << ", " << F_MDR << ", " << F_BULK << ", " << R << " || " << deltae1D << ", " << A << ", " << B << std::endl;

    //std::cout << gm->i << ", " << gm->j << ", " << (1.0/Vgeo[i])*Acon0[i]*delta_BULK*kappa*Ac << std::endl;

    // total force calculation
    (contactSide == 0) ? F0 = F_MDR + F_BULK : F1 = F_MDR + F_BULK;



    //std::cout << gm->i << ", " << gm->j << " | " << deltao << ", " << ddelta_bar[i] << ", " << R << ", " << psi[i] << ", " << psi_b << ", " << Ac << " | " << pij << ", " << wij << std::endl;

    // mean surface dipslacement calculation
     *Ac_offset = wij*Ac;

    // radius update scheme quantity calculation
    Vcaps[i] += (M_PI/3.0)*pow(delta,2.0)*(3.0*R - delta);
    
    const double Fntmp = wij*(F_MDR + F_BULK);
    const double fx = Fntmp*gm->nx[0];
    const double fy = Fntmp*gm->nx[1];
    const double fz = Fntmp*gm->nx[2];
    const double bx = -(Ro - deltao)*gm->nx[0];
    const double by = -(Ro - deltao)*gm->nx[1];
    const double bz = -(Ro - deltao)*gm->nx[2];
    const double eps_bar_contact = (1.0/(3*kappa*Velas[i]))*(fx*bx + fy*by + fz*bz);
    eps_bar[i] += eps_bar_contact;
    
      //if ( (i == 0 && j == 2 && gm->contact_type == 0) || (i == 2 && j == 0 && gm->contact_type == 0)) {
      //std::cout << i << ", " << j << ", " << gm->contact_type << " || " << delta << ", " << *delta_offset << ", " << (uintptr_t)(delta_offset) << " || " << deltao << ", " << *deltao_offset << ", " << (uintptr_t)(deltao_offset) << " || " << delta_MDR << ", " << *delta_MDR_offset << ", " << (uintptr_t)(delta_MDR_offset) << " || " << *Yflag_offset << ", " << (uintptr_t)(Yflag_offset) << " || " << R << std::endl;
      //std::cout << i << ", " << j << ", " << gm->contact_type << " || " << fx << ", " << fy << ", " << fz << " || " << bx << ", " << by << ", " << bz << ", " << Velas[i] << std::endl;
      //std::cout << i << ", " << j << ", " << gm->contact_type << " || " << eps_bar_contact << ", " << *eps_bar_offset << ", " << (uintptr_t)(eps_bar_offset) << " || " << wij << ", " << ddeltao << ", " << deltao << ", " << delta << ", " << *delta_offset << " || " << Ro << ", " << R << std::endl;
      //}

      //if () {
      //  std::cout << j << ", " << -Vo*(eps_bar_contact - *eps_bar_offset) - wij*M_PI*ddeltao*( 2.0*deltao*Ro - pow(deltao,2.0) + pow(R,2.0) - pow(Ro,2.0) ) << ", " << -Vo*(eps_bar_contact - *eps_bar_offset) << ", " << wij*M_PI*ddeltao*( 2.0*deltao*Ro - pow(deltao,2.0) + pow(R,2.0) - pow(Ro,2.0) ) << std::endl;
      //std::cout << i << ", " << j << ", " << gm->contact_type << " || " << eps_bar_contact << ", " << *eps_bar_offset << ", " << (uintptr_t)(eps_bar_offset) << " || " << wij << ", " << ddeltao << ", " << deltao << " || " << Ro << ", " << R << std::endl;
      //}

    if(delta_MDR == deltamax_MDR && *Yflag_offset > 0.0 && F_MDR > 0.0){
      const double Vo = (4.0/3.0)*M_PI*pow(Ro,3.0);
      dRnumerator[i] += -Vo*(eps_bar_contact - *eps_bar_offset) - wij*M_PI*ddeltao*( 2.0*deltao*Ro - pow(deltao,2.0) + pow(R,2.0) - pow(Ro,2.0) );

      dRdenominator[i] += wij*2.0*M_PI*R*(deltao + R - Ro);
    }
    *eps_bar_offset = eps_bar_contact;

    sigmaxx[i] += (1.0/Velas[i])*(fx*bx);
    sigmayy[i] += (1.0/Velas[i])*(fy*by);
    sigmazz[i] += (1.0/Velas[i])*(fz*bz);
    //std::cout << psi_b << ", " << psi[i] << ", " << A << ", " << B << ", " << pY << ", " << amax << " || " << deltao << ", " << delta << ", " << ddelta << ", " << *delta_offset << ", " << ddelta_bar[i] << " || " << delta_MDR << ", " << ddelta_MDR << ", " << *delta_MDR_offset << ", " << deltamax_MDR << " || " << delta_BULK << ", " << ddelta_BULK << ", " << *delta_BULK_offset << " || " << R << " || " << Ac << ", " << *Ac_offset << ", " << Acon0[i] << ", " << Acon1[i] << " || " << F_MDR << ", " << F_BULK << ", " << Vgeo[i] << std::endl;

    //std::cout << gm->i << ", " << gm->j << ", " << gm->radi << ", " << gm->radj << " | " << delta << ", " << F_MDR << " | " << deltae1D << ", " << A << ", " << B << std::endl;

  }

  // rigid flat placement scheme 
  *F_offset0 = F0;
  *F_offset1 = F1;

  gm->i = i_true;
  gm->j = j_true;
  gm->radi = radi_true;
  gm->radj = radj_true;

  double * penalty_offset = & history[penalty_offset_];
  const double pij = *penalty_offset;
  const double wij = std::max(1.0-pij,0.0);
  *penalty_offset = 0.0;

  //std::cout << gm->i << ", " << gm->j  << ", " << xi << ", " << xj << std::endl;

  // wall force magnifier
  double * deltao_offset = & history[deltao_offset_0];
  //const double wallForceMagnifer = std::exp(10.0*(*deltao_offset)/Rinitial[gm->i] - 9.0) + 1.0;
  const double wallForceMagnifer = 1.0;

  // assign final force
  //(gm->contact_type != PAIR) ? F = wij*F0*wallForceMagnifer : F = wij*(F0 + F1)/2;  // F = 2*wij*pow(1/F0 + 1/F1,-1);

  if (gm->contact_type != PAIR) {
    F = wij*F0*wallForceMagnifer;
  } else {
    F = wij*(F0 + F1)/2; 
    //if (R0 <  R1) {
    //  double b = pow(R1/R0,0.65);
    //  F = wij*F0; // wij*(F0 + F1/b)/2; //
    //  //F = wij*F0;
    //  //std::cout << "Number 1: " << F << ", " << F0 << ", " << F1 << " | " << R0 << ", " << R1 << ", " << b << std::endl; 
    //} else {
    //  double b = pow(R1/R0,0.65);
    //  F = wij*F1;  // wij*(F0/b + F1)/2; //
    //  //F = wij*F1;
    //  //std::cout << "Number 2: " << F << ", " << F0 << ", " << F1 << " | " << R0 << ", " << R1 << ", " << b << std::endl;
    //}
  }

  //std::cout << F << ", " << F0 << ", " << F1 << " | " << R0 << ", " << R1 << std::endl;

  // calculate damping force
  if (F > 0.0) {
    double Eeff;
    double Reff;
    if (gm->contact_type == PAIR) {
      Eeff = E/(2.0*(1.0-pow(nu,2.0)));
      Reff = pow((1/gm->radi + 1/gm->radj),-1);
    } else {
      Eeff = E/(1.0-pow(nu,2.0));
      Reff = gm->radi;
    }
    const double kn = Eeff*Reff;
    const double beta = -log(CoR)/sqrt(pow(log(CoR),2.0) + M_PI*M_PI);
    const double damp_prefactor = beta*sqrt(gm->meff*kn);
    const double F_DAMP = -damp_prefactor*(gm->vnnr);

    //std:: cout << gm->contact_type << ", " << Eeff << " , " << Reff << ", " << gm->radi << ", " << gm->radj << " || " << kn << ", " << beta << ", " << gm->meff << " || " << F_DAMP << ", " << F << std::endl;
    F += wij*F_DAMP;
  }

  //double **x = atom->x;
  //const double xi = x[gm->i][0];
  //const double xj = x[gm->j][0];
  //const double del = 20.0 - abs(xi-xj);
  
  //if (gm->i == 0 && gm->j == 1) {
  //  CSVWriter csvWriter("/Users/willzunker/lammps/sims/twoParticleDifferingRadii/contact_radius.csv");
  //  std::stringstream rowDataStream;
  //  rowDataStream << std::scientific << std::setprecision(6); // Set the format and precision
  //  rowDataStream << a0 << ", " << a1;
  //  std::string rowData = rowDataStream.str();
  //  csvWriter.writeRow(rowData);
  //}

  //if (gm->i == 0 && gm->j == 2) {
  //  CSVWriter csvWriter("/Users/willzunker/lammps/sims/compressionSleeve/pairContactsBotCen.csv");
  //  std::stringstream rowDataStream;
  //  rowDataStream << std::scientific << std::setprecision(4); // Set the format and precision
  //  rowDataStream << del << ", " << F;
  //  std::string rowData = rowDataStream.str();
  //  csvWriter.writeRow(rowData);
  //}

  return F;
}

//std::cout << sidata.i << ", " << sidata.j << ", " << R << ", " << deltan << ", " << deltao << ", " << dRsums_i[0] << ", " << dRsums_i[1] << ", " << numQuant << std::endl;
//std::cout << gm->i << ", " << gm->j << " | " << gm->nx[0] << ", " << gm->nx[1] << ", " << gm->nx[2] << std::endl;
//std::cout << "F is: " << F << std::endl;
//std::cout << "Ro from particle history is: " << Ro[gm->i] << std::endl;
//std::cout << "MDR contact model has been entered." << std::endl;
// std::cout << "F is: " << F << std::endl;
//std::cout << "gamma > 0.0: " << F_MDR << ", " << gm->i << ", " << gm->j << std::endl;
//std::cout << "deltae1D <= 0.0: " << F_MDR << ", " << gm->i << ", " << gm->j << std::endl;
//std::cout << "F_MDR should be > 0: " << F_MDR << ", " << gm->i << ", " << gm->j << std::endl;
//std::cout << gm->i << ", " << gm->j << " || " << delta << ", " << delta_MDR << ", " << deltamax_MDR << ", " << deltae1D << " || " << A << ", " << B << std::endl;

//std::cout << i_true << ", " << j_true << std::endl; 

      //std::cout << history_index << ", " << history[0] << ", " << history[1] << ", " << history[2] << std::endl;

      // initialize all history variables
      //double delta_offset; 
      //double deltao_offset;
      //double delta_MDR_offset;   
      //double delta_BULK_offset; 
      //double deltamax_MDR_offset; 
      //double Yflag; 
      //double deltaY_offset; 
      //double Ac_offset; 
      //double aAdh_offset; 
      //double deltap_offset; 
      //double cA_offset; 
      //double eps_bar_offset;
      //double wall_contact_flag_offset;