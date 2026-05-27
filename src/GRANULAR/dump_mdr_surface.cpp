/* ----------------------------------------------------------------------
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
   Contributing authors:
   Dalil Ashong (UC Berkeley), William Zunker (MIT), Ken Kamrin (UC Berkeley)
------------------------------------------------------------------------- */

#include "dump_mdr_surface.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "fix_granular_mdr.h"
#include "fix_neigh_history.h"
#include "fix_wall_gran.h"
#include "fix_wall_gran_region.h"
#include "force.h"
#include "gran_sub_mod_normal.h"
#include "granular_model.h"
#include "mdr_reconstruct.h"
#include "mdr_surface_recon.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "pair_granular.h"
#include "region.h"
#include "update.h"
#include "utils.h"

#include <cmath>
#include <cstdint>
#include <cstring>

using namespace LAMMPS_NS;
using namespace Granular_NS;
using namespace Granular_MDR_NS;

/* ---------------------------------------------------------------------- */

DumpMDRSurface::DumpMDRSurface(LAMMPS *lmp, int narg, char **arg) : Dump(lmp, narg, arg)
{
  if (narg < 5) error->all(FLERR, "Illegal dump mdr/surface command");

  ntheta = 60;
  nphi = 30;
  csv_flag = 0;
  csvfile = nullptr;
  pair = nullptr;
  mdr_model = nullptr;
  fix_hist = nullptr;

  // one mesh file per snapshot is required (filename must contain '*')
  if (!multifile)
    error->all(FLERR, "Dump mdr/surface requires one snapshot per file (use '*' in filename)");

  int iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "ntheta") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal dump mdr/surface command");
      ntheta = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "nphi") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal dump mdr/surface command");
      nphi = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "csv") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal dump mdr/surface command");
      csv_flag = 1;
      csvfile = utils::strdup(arg[iarg + 1]);
      iarg += 2;
    } else {
      error->all(FLERR, "Unknown dump mdr/surface keyword: {}", arg[iarg]);
    }
  }

  if (ntheta < 3 || nphi < 2)
    error->all(FLERR, "Dump mdr/surface requires ntheta >= 3 and nphi >= 2");

  // derive pvd filename from the dump filename (strip '*', swap suffix)
  std::string base = filename;
  auto star = base.find('*');
  if (star != std::string::npos) base.erase(star, 1);
  auto dot = base.rfind('.');
  if (dot != std::string::npos) base.erase(dot);
  pvdname = base + ".pvd";
}

/* ---------------------------------------------------------------------- */

DumpMDRSurface::~DumpMDRSurface()
{
  delete[] csvfile;
}

/* ---------------------------------------------------------------------- */

void DumpMDRSurface::find_mdr_model()
{
  pair = dynamic_cast<PairGranular *>(force->pair_match("granular", 1));
  if (!pair) error->all(FLERR, "Dump mdr/surface requires pair_style granular with the MDR model");

  for (int i = 0; i < pair->nmodels; i++) {
    GranularModel *gm = pair->models_list[i];
    if (gm->normal_model && gm->normal_model->name == "mdr") {
      mdr_model = dynamic_cast<GranSubModNormalMDR *>(gm->normal_model);
      break;
    }
  }
  if (!mdr_model) error->all(FLERR, "Dump mdr/surface requires an MDR normal model");

  fix_hist = dynamic_cast<FixNeighHistory *>(modify->get_fix_by_id("NEIGH_HISTORY_GRANULAR"));
  if (!fix_hist) error->all(FLERR, "Dump mdr/surface cannot find granular contact history");
}

/* ---------------------------------------------------------------------- */

void DumpMDRSurface::init_style()
{
  if (binary) error->all(FLERR, "Dump mdr/surface does not support binary files");
  find_mdr_model();

  int dim, cols;
  index_Ro = atom->find_custom("Ro", dim, cols);
  index_psi = atom->find_custom("psi", dim, cols);
  index_sigmaxx = atom->find_custom("sigmaxx", dim, cols);
  index_sigmayy = atom->find_custom("sigmayy", dim, cols);
  index_sigmazz = atom->find_custom("sigmazz", dim, cols);
}

/* ----------------------------------------------------------------------
   per-step filename with '*' replaced by the timestep
------------------------------------------------------------------------- */

std::string DumpMDRSurface::this_step_name(const char *templatename)
{
  char *fc = utils::strdup(utils::star_subst(templatename, update->ntimestep, padflag));
  std::string s(fc);
  delete[] fc;
  return s;
}

/* ---------------------------------------------------------------------- */

void DumpMDRSurface::write()
{
  double **x = atom->x;
  double *radius = atom->radius;
  tagint *tag = atom->tag;
  int *mask = atom->mask;
  const int nlocal = atom->nlocal;

  double *psi = (index_psi >= 0) ? atom->dvector[index_psi] : nullptr;
  double *sxx = (index_sigmaxx >= 0) ? atom->dvector[index_sigmaxx] : nullptr;
  double *syy = (index_sigmayy >= 0) ? atom->dvector[index_sigmayy] : nullptr;
  double *szz = (index_sigmazz >= 0) ? atom->dvector[index_sigmazz] : nullptr;
  double *Rinit = (index_Ro >= 0) ? atom->dvector[index_Ro] : nullptr;
  std::vector<double> csv_rows;    // 21 cols/row: tag nx ny nz Ro R delta deltae A B a_na deltamax E nu Y deltaR amax x y z step

  const double E = mdr_model->get_emod();
  const double nu = mdr_model->get_poiss();
  const double Y = mdr_model->Y;
  const double Eeff = E / (1.0 - nu * nu);
  const double Eeffinv = 1.0 / Eeff;
  const double G = E / (2.0 * (1.0 + nu));

  const int stride = pair->get_size_history();
  NeighList *list = pair->list;
  const int inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  double **firsthistory = fix_hist->firstvalue;

  // ---- collect caps per owned atom (pair contacts), dedup by partner tag ----
  std::vector<std::vector<MDRCap>> caps(nlocal);
  std::vector<std::vector<tagint>> cappart(nlocal);

  // core: build one cap on side atom s from history block h (side index 0/1),
  // with outward normal n pointing from the particle centre toward the contact
  auto build_cap = [&](int s, double *h, int side, double nx, double ny, double nz) {
    const double yflag = h[YFLAG_0 + side];
    const double delta = h[DELTA_0 + side];
    if (delta <= 1e-15) return;    // no current overlap
    const double delta_MDR = h[DELTA_MDR_0 + side];
    const double deltamax_MDR = h[DELTAMAX_MDR_0 + side];
    const double cA = h[CA_0 + side];
    const double deltap = h[DELTAP_0 + side];
    const double deltamax_app = h[DELTA_MAX];
    const double R = radius[s];

    const double pY = mdr_pressure_yield(Y, deltamax_MDR, R);
    const MDRCapGeometry g =
        mdr_cap_geometry(yflag, R, delta_MDR, deltamax_MDR, cA, pY, Eeff, Eeffinv, G, nu);
    const double deltae = deltamax_MDR - deltap;

    MDRCap c;
    c.n[0] = nx; c.n[1] = ny; c.n[2] = nz;
    c.profile = mdr_cap_profile(yflag, R, delta, deltae, g.A, g.B, g.a_na, deltamax_app, g.deltaR,
                                g.amax);
    caps[s].push_back(c);

    // remember the deepest plastic imprint for this contact direction so it can
    // be redrawn (elastically sprung back) after the contact separates
    if (yflag != 0.0) {
      std::vector<PersistCap> &v = persist[tag[s]];
      int found = -1;
      for (size_t k = 0; k < v.size(); k++)
        if (v[k].n[0] * nx + v[k].n[1] * ny + v[k].n[2] * nz > 0.951) { found = (int) k; break; }
      if (found < 0 || deltamax_app > v[found].deltamax) {
        PersistCap pc;
        pc.n[0] = nx; pc.n[1] = ny; pc.n[2] = nz;
        pc.R = R; pc.delta = delta; pc.deltae = deltae; pc.A = g.A; pc.B = g.B;
        pc.a_na = g.a_na; pc.deltamax = deltamax_app; pc.deltaR = g.deltaR; pc.amax = g.amax;
        if (found < 0) v.push_back(pc);
        else v[found] = pc;
      }
    }

    if (csv_flag) {
      const double Ro = Rinit ? Rinit[s] : R;
      const double row[21] = {(double) tag[s], nx, ny, nz, Ro, R, delta, deltae,
                              g.A, g.B, g.a_na, deltamax_app, E, nu, Y, g.deltaR, g.amax,
                              x[s][0], x[s][1], x[s][2], (double) update->ntimestep};
      for (double rv : row) csv_rows.push_back(rv);
    }
  };

  // pair contacts: normal points from the side atom toward its partner
  auto add_pair_cap = [&](int s, int o, double *h, int side) {
    for (auto pt : cappart[s])
      if (pt == tag[o]) return;    // already have this contact
    double dx = x[o][0] - x[s][0], dy = x[o][1] - x[s][1], dz = x[o][2] - x[s][2];
    double r = std::sqrt(dx * dx + dy * dy + dz * dz);
    if (r < 1e-15) return;
    double rinv = 1.0 / r;
    cappart[s].push_back(tag[o]);
    build_cap(s, h, side, dx * rinv, dy * rinv, dz * rinv);
  };

  for (int ii = 0; ii < inum; ii++) {
    const int i = ilist[ii];
    double *allh = firsthistory[i];
    const int jnum = numneigh[i];
    int *jlist = firstneigh[i];
    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj] & NEIGHMASK;
      double radsum = radius[i] + radius[j];
      double dx = x[j][0] - x[i][0], dy = x[j][1] - x[i][1], dz = x[j][2] - x[i][2];
      if (dx * dx + dy * dy + dz * dz >= radsum * radsum) continue;    // not touching
      double *h = &allh[stride * jj];
      if (i < nlocal && (mask[i] & groupbit)) add_pair_cap(i, j, h, (tag[i] > tag[j]) ? 0 : 1);
      if (j < nlocal && (mask[j] & groupbit)) add_pair_cap(j, i, h, (tag[j] > tag[i]) ? 0 : 1);
    }
  }

  // wall/gran/region contacts: normal points from the particle toward the wall.
  // Re-query each region for geometry (positions unchanged since the force eval)
  // and match to the persisted per-contact history by wall index.
  for (auto *wf : modify->get_fix_by_style("wall/gran/region")) {
    auto *fwr = dynamic_cast<FixWallGranRegion *>(wf);
    if (!fwr || !fwr->model || !fwr->model->normal_model) continue;
    if (fwr->model->normal_model->name != "mdr") continue;
    Region *reg = fwr->region;
    reg->prematch();
    for (int i = 0; i < nlocal; i++) {
      if (!(mask[i] & groupbit)) continue;
      const int nc = reg->surface(x[i][0], x[i][1], x[i][2], radius[i]);
      for (int ic = 0; ic < nc; ic++) {
        const int iwall = reg->contact[ic].iwall;
        int slot = -1;
        for (int c = 0; c < fwr->ncontact[i]; c++)
          if (fwr->walls[i][c] == iwall) { slot = c; break; }
        if (slot < 0) continue;
        const double rr = reg->contact[ic].r;
        if (rr < 1e-15) continue;
        const double rinv = -1.0 / rr;    // contact.del points wall->particle; flip to face wall
        build_cap(i, fwr->history_many[i][slot], 0, reg->contact[ic].delx * rinv,
                  reg->contact[ic].dely * rinv, reg->contact[ic].delz * rinv);
      }
    }
  }

  // ---- persistent residual imprints ----
  // For every plastic contact ever made, add its deepest cap sprung back outward
  // by the elastic recovery deltae. The min-over-caps then picks the deeper of
  // the live cap (when still loaded, delta > deltap) and this residual (once
  // unloaded past the residual, or fully separated). This makes recovery stop
  // continuously at the elastic spring-back instead of rounding to a sphere and
  // then snapping back when the contact leaves the neighbor list.
  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    auto it = persist.find(tag[i]);
    if (it == persist.end()) continue;
    for (const PersistCap &pc : it->second) {
      MDRCapProfile deep = mdr_cap_profile(1.0, pc.R, pc.delta, pc.deltae, pc.A, pc.B, pc.a_na,
                                           pc.deltamax, pc.deltaR, pc.amax);
      MDRCap c;
      c.n[0] = pc.n[0]; c.n[1] = pc.n[1]; c.n[2] = pc.n[2];
      c.profile = mdr_lift_profile(deep, pc.deltae);    // elastic spring-back
      if (!c.profile.alpha.empty()) caps[i].push_back(c);
    }
  }

  // ---- build meshes for owned group atoms; serialize for this proc ----
  std::vector<double> pxyz;        // 3 per vertex
  std::vector<int> pconn;          // 3 per triangle (proc-local vertex idx)
  std::vector<double> cell_id, cell_R, cell_psi, cell_vm;    // per triangle

  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    const double R = radius[i];
    MDRMesh m = mdr_build_particle_mesh(x[i], R, caps[i], ntheta, nphi);
    const int voff = (int) (pxyz.size() / 3);
    for (double v : m.xyz) pxyz.push_back(v);
    const int ntri = (int) (m.tris.size() / 3);
    double vm = 0.0;
    if (sxx && syy && szz) {
      double a = sxx[i] - syy[i], b = syy[i] - szz[i], c = szz[i] - sxx[i];
      vm = std::sqrt(0.5 * (a * a + b * b + c * c));
    }
    const double pv = psi ? psi[i] : 0.0;
    for (int t = 0; t < ntri; t++) {
      pconn.push_back(m.tris[3 * t + 0] + voff);
      pconn.push_back(m.tris[3 * t + 1] + voff);
      pconn.push_back(m.tris[3 * t + 2] + voff);
      cell_id.push_back((double) tag[i]);
      cell_R.push_back(R);
      cell_psi.push_back(pv);
      cell_vm.push_back(vm);
    }
  }

  // ---- gather to proc 0 ----
  int my_nv = (int) (pxyz.size() / 3);
  int my_nt = (int) (pconn.size() / 3);
  std::vector<int> nv_all(nprocs), nt_all(nprocs);
  MPI_Gather(&my_nv, 1, MPI_INT, nv_all.data(), 1, MPI_INT, 0, world);
  MPI_Gather(&my_nt, 1, MPI_INT, nt_all.data(), 1, MPI_INT, 0, world);

  std::vector<int> v_disp(nprocs, 0), t_disp(nprocs, 0);
  std::vector<int> v3_cnt(nprocs, 0), v3_disp(nprocs, 0), t3_cnt(nprocs, 0), t3_disp(nprocs, 0);
  int tot_nv = 0, tot_nt = 0;
  if (me == 0) {
    for (int p = 0; p < nprocs; p++) {
      v_disp[p] = tot_nv; t_disp[p] = tot_nt;
      tot_nv += nv_all[p]; tot_nt += nt_all[p];
    }
    for (int p = 0; p < nprocs; p++) {
      v3_cnt[p] = 3 * nv_all[p]; v3_disp[p] = 3 * v_disp[p];
      t3_cnt[p] = 3 * nt_all[p]; t3_disp[p] = 3 * t_disp[p];
    }
  }

  std::vector<double> all_xyz, all_id, all_R, all_psi, all_vm;
  std::vector<int> all_conn;
  if (me == 0) {
    all_xyz.resize(3 * tot_nv);
    all_conn.resize(3 * tot_nt);
    all_id.resize(tot_nt); all_R.resize(tot_nt); all_psi.resize(tot_nt); all_vm.resize(tot_nt);
  }
  MPI_Gatherv(pxyz.data(), 3 * my_nv, MPI_DOUBLE, all_xyz.data(), v3_cnt.data(), v3_disp.data(),
              MPI_DOUBLE, 0, world);
  MPI_Gatherv(pconn.data(), 3 * my_nt, MPI_INT, all_conn.data(), t3_cnt.data(), t3_disp.data(),
              MPI_INT, 0, world);
  MPI_Gatherv(cell_id.data(), my_nt, MPI_DOUBLE, all_id.data(), nt_all.data(), t_disp.data(),
              MPI_DOUBLE, 0, world);
  MPI_Gatherv(cell_R.data(), my_nt, MPI_DOUBLE, all_R.data(), nt_all.data(), t_disp.data(),
              MPI_DOUBLE, 0, world);
  MPI_Gatherv(cell_psi.data(), my_nt, MPI_DOUBLE, all_psi.data(), nt_all.data(), t_disp.data(),
              MPI_DOUBLE, 0, world);
  MPI_Gatherv(cell_vm.data(), my_nt, MPI_DOUBLE, all_vm.data(), nt_all.data(), t_disp.data(),
              MPI_DOUBLE, 0, world);

  // offset each proc's connectivity by its global vertex base
  if (me == 0) {
    for (int p = 0; p < nprocs; p++)
      for (int k = 0; k < 3 * nt_all[p]; k++) all_conn[3 * t_disp[p] + k] += v_disp[p];
  }

  // ---- write VTP on proc 0 (binary, appended raw) ----
  if (me == 0) {
    std::string fname = this_step_name(filename);
    FILE *f = fopen(fname.c_str(), "wb");
    if (!f) error->one(FLERR, "Cannot open dump mdr/surface file {}", fname);

    // pack little-endian typed buffers (native = little-endian on x86/ARM)
    std::vector<float> b_pts(3 * tot_nv);
    for (int k = 0; k < 3 * tot_nv; k++) b_pts[k] = (float) all_xyz[k];
    std::vector<int32_t> b_conn(3 * tot_nt);
    for (int k = 0; k < 3 * tot_nt; k++) b_conn[k] = (int32_t) all_conn[k];
    std::vector<int32_t> b_off(tot_nt);
    for (int t = 0; t < tot_nt; t++) b_off[t] = 3 * (t + 1);
    std::vector<int32_t> b_id(tot_nt);
    for (int t = 0; t < tot_nt; t++) b_id[t] = (int32_t) all_id[t];
    std::vector<float> b_rad(tot_nt), b_psi(tot_nt), b_vm(tot_nt);
    for (int t = 0; t < tot_nt; t++) {
      b_rad[t] = (float) all_R[t];
      b_psi[t] = (float) all_psi[t];
      b_vm[t] = (float) all_vm[t];
    }

    // appended blob: each array prefixed by a UInt64 byte count; offsets are
    // measured from the first byte after the leading underscore
    std::string blob;
    auto put = [&](const void *d, size_t nbytes) -> uint64_t {
      const uint64_t off = blob.size();
      const uint64_t n = nbytes;
      blob.append(reinterpret_cast<const char *>(&n), sizeof(uint64_t));
      blob.append(reinterpret_cast<const char *>(d), nbytes);
      return off;
    };
    const uint64_t o_pts = put(b_pts.data(), b_pts.size() * sizeof(float));
    const uint64_t o_conn = put(b_conn.data(), b_conn.size() * sizeof(int32_t));
    const uint64_t o_off = put(b_off.data(), b_off.size() * sizeof(int32_t));
    const uint64_t o_id = put(b_id.data(), b_id.size() * sizeof(int32_t));
    const uint64_t o_rad = put(b_rad.data(), b_rad.size() * sizeof(float));
    uint64_t o_psi = 0, o_vm = 0;
    if (psi) o_psi = put(b_psi.data(), b_psi.size() * sizeof(float));
    if (sxx) o_vm = put(b_vm.data(), b_vm.size() * sizeof(float));

    fprintf(f, "<?xml version=\"1.0\"?>\n");
    fprintf(f, "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\" "
               "header_type=\"UInt64\">\n <PolyData>\n");
    fprintf(f,
            "  <Piece NumberOfPoints=\"%d\" NumberOfVerts=\"0\" NumberOfLines=\"0\" "
            "NumberOfStrips=\"0\" NumberOfPolys=\"%d\">\n",
            tot_nv, tot_nt);
    fprintf(f,
            "   <Points>\n    <DataArray type=\"Float32\" NumberOfComponents=\"3\" "
            "format=\"appended\" offset=\"%llu\"/>\n   </Points>\n",
            (unsigned long long) o_pts);
    fprintf(f, "   <Polys>\n");
    fprintf(f, "    <DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" "
               "offset=\"%llu\"/>\n", (unsigned long long) o_conn);
    fprintf(f, "    <DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" "
               "offset=\"%llu\"/>\n   </Polys>\n", (unsigned long long) o_off);
    fprintf(f, "   <CellData Scalars=\"particle_id\">\n");
    fprintf(f, "    <DataArray type=\"Int32\" Name=\"particle_id\" format=\"appended\" "
               "offset=\"%llu\"/>\n", (unsigned long long) o_id);
    fprintf(f, "    <DataArray type=\"Float32\" Name=\"radius\" format=\"appended\" "
               "offset=\"%llu\"/>\n", (unsigned long long) o_rad);
    if (psi)
      fprintf(f, "    <DataArray type=\"Float32\" Name=\"psi\" format=\"appended\" "
                 "offset=\"%llu\"/>\n", (unsigned long long) o_psi);
    if (sxx)
      fprintf(f, "    <DataArray type=\"Float32\" Name=\"von_mises\" format=\"appended\" "
                 "offset=\"%llu\"/>\n", (unsigned long long) o_vm);
    fprintf(f, "   </CellData>\n  </Piece>\n </PolyData>\n");

    fprintf(f, " <AppendedData encoding=\"raw\">\n  _");
    fwrite(blob.data(), 1, blob.size(), f);
    fprintf(f, "\n </AppendedData>\n</VTKFile>\n");
    fclose(f);

    // record + rewrite pvd time series (basename for portability)
    std::string bn = fname;
    auto slash = bn.find_last_of("/\\");
    if (slash != std::string::npos) bn = bn.substr(slash + 1);
    pvd_entries.emplace_back((double) update->ntimestep, bn);
    write_pvd();
  }

  // ---- optional raw per-contact CSV (debug / MATLAB-compatible) ----
  if (csv_flag) {
    int my_n = (int) csv_rows.size();    // multiple of 21
    std::vector<int> cnt(nprocs), disp(nprocs, 0);
    MPI_Gather(&my_n, 1, MPI_INT, cnt.data(), 1, MPI_INT, 0, world);
    int tot = 0;
    std::vector<double> allrows;
    if (me == 0) {
      for (int p = 0; p < nprocs; p++) { disp[p] = tot; tot += cnt[p]; }
      allrows.resize(tot);
    }
    MPI_Gatherv(csv_rows.data(), my_n, MPI_DOUBLE, allrows.data(), cnt.data(), disp.data(),
                MPI_DOUBLE, 0, world);
    if (me == 0) {
      std::string cname = this_step_name(csvfile);
      FILE *cf = fopen(cname.c_str(), "w");
      if (cf) {
        fprintf(cf,
                "tag,nx,ny,nz,Ro,R,delta,deltae,A,B,a_na,deltamax,E,nu,Y,deltaR,amax,x,y,z,step\n");
        const int nrow = tot / 21;
        for (int r = 0; r < nrow; r++) {
          for (int k = 0; k < 21; k++) fprintf(cf, k ? ",%.10g" : "%.10g", allrows[r * 21 + k]);
          fprintf(cf, "\n");
        }
        fclose(cf);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void DumpMDRSurface::write_pvd()
{
  FILE *f = fopen(pvdname.c_str(), "w");
  if (!f) return;
  fprintf(f, "<?xml version=\"1.0\"?>\n");
  fprintf(f, "<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\">\n");
  fprintf(f, " <Collection>\n");
  for (auto &e : pvd_entries)
    fprintf(f, "  <DataSet timestep=\"%g\" group=\"\" part=\"0\" file=\"%s\"/>\n", e.first,
            e.second.c_str());
  fprintf(f, " </Collection>\n</VTKFile>\n");
  fclose(f);
}
