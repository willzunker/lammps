/***************************************************************************
                               base_amoeba.cpp
                             -------------------
                            Trung Dac Nguyen (Northwestern)

  Base class for pair styles needing per-particle data for position,
  charge, and type.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                :
    email                : trung.nguyen@northwestern.edu
 ***************************************************************************/

#include "lal_base_amoeba.h"
namespace LAMMPS_AL {
#define BaseAmoebaT BaseAmoeba<numtyp, acctyp>

extern Device<PRECISION,ACC_PRECISION> global_device;

template <class numtyp, class acctyp>
BaseAmoebaT::BaseAmoeba() : _compiled(false), _max_bytes(0) {
  device=&global_device;
  ans=new Answer<numtyp,acctyp>();
  nbor=new Neighbor();
  pair_program=nullptr;
  ucl_device=nullptr;
  #if defined(LAL_OCL_EV_JIT)
  pair_program_noev=nullptr;
  #endif
}

template <class numtyp, class acctyp>
BaseAmoebaT::~BaseAmoeba() {
  delete ans;
  delete nbor;
  k_polar.clear();
  k_special15.clear();
  if (pair_program) delete pair_program;
}

template <class numtyp, class acctyp>
int BaseAmoebaT::bytes_per_atom_atomic(const int max_nbors) const {
  return device->atom.bytes_per_atom()+ans->bytes_per_atom()+
         nbor->bytes_per_atom(max_nbors);
}

template <class numtyp, class acctyp>
int BaseAmoebaT::init_atomic(const int nlocal, const int nall,
                             const int max_nbors, const int maxspecial,
                             const int maxspecial15,
                             const double cell_size, const double gpu_split,
                             FILE *_screen, const void *pair_program,
                             const char *k_name) {
  screen=_screen;

  int gpu_nbor=0;
  if (device->gpu_mode()==Device<numtyp,acctyp>::GPU_NEIGH)
    gpu_nbor=1;
  else if (device->gpu_mode()==Device<numtyp,acctyp>::GPU_HYB_NEIGH)
    gpu_nbor=2;

  int _gpu_host=0;
  int host_nlocal=hd_balancer.first_host_count(nlocal,gpu_split,gpu_nbor);
  if (host_nlocal>0)
    _gpu_host=1;

  _threads_per_atom=device->threads_per_charge();

  bool charge = true;
  bool rot = false;
  bool vel = false;
  _extra_fields = 24; // round up to accomodate quadruples of numtyp values
                      // rpole 13; uind 3; uinp 3; amtype, amgroup
  int success=device->init(*ans,charge,rot,nlocal,nall,maxspecial,vel,_extra_fields);
  if (success!=0)
    return success;

  if (ucl_device!=device->gpu) _compiled=false;

  ucl_device=device->gpu;
  atom=&device->atom;

  _block_size=device->pair_block_size();
  _block_bio_size=device->block_bio_pair();
  compile_kernels(*ucl_device,pair_program,k_name);

  if (_threads_per_atom>1 && gpu_nbor==0) {
    nbor->packing(true);
    _nbor_data=&(nbor->dev_packed);
  } else
    _nbor_data=&(nbor->dev_nbor);

  success = device->init_nbor(nbor,nlocal,host_nlocal,nall,maxspecial,_gpu_host,
                  max_nbors,cell_size,false,_threads_per_atom);
  if (success!=0)
    return success;

  // Initialize host-device load balancer
  hd_balancer.init(device,gpu_nbor,gpu_split);

  // Initialize timers for the selected GPU
  time_pair.init(*ucl_device);
  time_pair.zero();

  pos_tex.bind_float(atom->x,4);
  q_tex.bind_float(atom->q,1);

  _max_an_bytes=ans->gpu_bytes()+nbor->gpu_bytes();
  _maxspecial=maxspecial;
  _maxspecial15=maxspecial15;

  // allocate per-atom array tep 

  int ef_nall=nall;
  if (ef_nall==0)
    ef_nall=2000;

  _max_tep_size=static_cast<int>(static_cast<double>(ef_nall)*1.10);
  _tep.alloc(_max_tep_size*4,*(this->ucl_device),UCL_READ_WRITE,UCL_READ_WRITE);
  dev_nspecial15.alloc(nall,*(this->ucl_device),UCL_READ_ONLY);
  dev_special15.alloc(_maxspecial15*nall,*(this->ucl_device),UCL_READ_ONLY);
  dev_special15_t.alloc(nall*_maxspecial15,*(this->ucl_device),UCL_READ_ONLY);

  return success;
}

template <class numtyp, class acctyp>
void BaseAmoebaT::estimate_gpu_overhead(const int add_kernels) {
  device->estimate_gpu_overhead(1+add_kernels,_gpu_overhead,_driver_overhead);
}

template <class numtyp, class acctyp>
void BaseAmoebaT::clear_atomic() {
  // Output any timing information
  acc_timers();
  double avg_split=hd_balancer.all_avg_split();
  _gpu_overhead*=hd_balancer.timestep();
  _driver_overhead*=hd_balancer.timestep();
  device->output_times(time_pair,*ans,*nbor,avg_split,_max_bytes+_max_an_bytes,
                       _gpu_overhead,_driver_overhead,_threads_per_atom,screen);

  time_pair.clear();
  hd_balancer.clear();

  nbor->clear();
  ans->clear();

  _tep.clear();
  dev_nspecial15.clear();
  dev_special15.clear();
  dev_special15_t.clear();
}

// ---------------------------------------------------------------------------
// Copy neighbor list from host
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
int * BaseAmoebaT::reset_nbors(const int nall, const int inum, int *ilist,
                                   int *numj, int **firstneigh, bool &success) {
  success=true;

  int mn=nbor->max_nbor_loop(inum,numj,ilist);
  resize_atom(inum,nall,success);
  resize_local(inum,mn,success);
  if (!success)
    return nullptr;

  nbor->get_host(inum,ilist,numj,firstneigh,block_size());

  double bytes=ans->gpu_bytes()+nbor->gpu_bytes();
  if (bytes>_max_an_bytes)
    _max_an_bytes=bytes;

  return ilist;
}

// ---------------------------------------------------------------------------
// Build neighbor list on device
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
inline void BaseAmoebaT::build_nbor_list(const int inum, const int host_inum,
                                         const int nall, double **host_x,
                                         int *host_type, double *sublo,
                                         double *subhi, tagint *tag,
                                         int **nspecial, tagint **special,
                                         int *nspecial15, tagint **special15,
                                         bool &success) {
  success=true;
  resize_atom(inum,nall,success);
  resize_local(inum,host_inum,nbor->max_nbors(),success);
  if (!success)
    return;
  atom->cast_copy_x(host_x,host_type);

  int mn;
  nbor->build_nbor_list(host_x, inum, host_inum, nall, *atom, sublo, subhi,
                        tag, nspecial, special, success, mn, ans->error_flag);

  // add one-five neighbors

  if (_maxspecial15>0) {
    UCL_H_Vec<int> view_nspecial15;
    UCL_H_Vec<tagint> view_special15;
    view_nspecial15.view(nspecial15,nall,*ucl_device);
    view_special15.view(special15[0],nall*_maxspecial15,*ucl_device);
    ucl_copy(dev_nspecial15,view_nspecial15,nall,false);
    ucl_copy(dev_special15_t,view_special15,_maxspecial15*nall,false);
    nbor->transpose(dev_special15, dev_special15_t, _maxspecial15, nall);

    add_onefive_neighbors();
  }

  double bytes=ans->gpu_bytes()+nbor->gpu_bytes();
  if (bytes>_max_an_bytes)
    _max_an_bytes=bytes;
}

// ---------------------------------------------------------------------------
// Copy nbor list from host if necessary and then calculate forces, virials,..
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
void BaseAmoebaT::compute(const int f_ago, const int inum_full, const int nall,
                          double **host_x, int *host_type, int *host_amtype,
                          int *host_amgroup, double **host_rpole,
                          double **host_uind, double **host_uinp,
                          int *ilist, int *numj, int **firstneigh,
                          const bool eflag_in, const bool vflag_in,
                          const bool eatom, const bool vatom,
                          int &host_start, const double cpu_time,
                          bool &success, double *host_q, const int nlocal,
                          double *boxlo, double *prd, void **tep_ptr) {
  acc_timers();
  int eflag, vflag;
  if (eatom) eflag=2;
  else if (eflag_in) eflag=1;
  else eflag=0;
  if (vatom) vflag=2;
  else if (vflag_in) vflag=1;
  else vflag=0;

  #ifdef LAL_NO_BLOCK_REDUCE
  if (eflag) eflag=2;
  if (vflag) vflag=2;
  #endif

  set_kernel(eflag,vflag);

  // ------------------- Resize _tep array ------------------------

  if (nall>_max_tep_size) {
    _max_tep_size=static_cast<int>(static_cast<double>(nall)*1.10);
    _tep.resize(_max_tep_size*4);

    dev_nspecial15.clear();
    dev_special15.clear();
    dev_special15_t.clear();
    dev_nspecial15.alloc(nall,*(this->ucl_device),UCL_READ_ONLY);
    dev_special15.alloc(_maxspecial15*nall,*(this->ucl_device),UCL_READ_ONLY);
    dev_special15_t.alloc(nall*_maxspecial15,*(this->ucl_device),UCL_READ_ONLY);   
  }

  *tep_ptr=_tep.host.begin();

  if (inum_full==0) {
    host_start=0;
    // Make sure textures are correct if realloc by a different hybrid style
    resize_atom(0,nall,success);
    zero_timers();
    return;
  }

  int ago=hd_balancer.ago_first(f_ago);
  int inum=hd_balancer.balance(ago,inum_full,cpu_time);
  ans->inum(inum);
  host_start=inum;

  if (ago==0) {
    reset_nbors(nall, inum, ilist, numj, firstneigh, success);
    if (!success)
      return;
  }

  // packing host arrays into host_extra

  atom->cast_x_data(host_x,host_type);
  atom->cast_q_data(host_q);
  cast_extra_data(host_amtype, host_amgroup, host_rpole, host_uind, host_uinp);
  hd_balancer.start_timer();
  atom->add_x_data(host_x,host_type);
  atom->add_q_data();
  atom->add_extra_data();

  device->precompute(f_ago,nlocal,nall,host_x,host_type,success,host_q,
                     boxlo, prd);

  const int red_blocks=loop(eflag,vflag);
  ans->copy_answers(eflag_in,vflag_in,eatom,vatom,ilist,red_blocks);
  device->add_ans_object(ans);
  hd_balancer.stop_timer();
}

// ---------------------------------------------------------------------------
// Reneighbor on GPU if necessary and then compute forces, virials, energies
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
int** BaseAmoebaT::compute(const int ago, const int inum_full, const int nall,
                           double **host_x, int *host_type, int *host_amtype,
                           int *host_amgroup, double **host_rpole,
                           double **host_uind, double **host_uinp,
                           double *sublo, double *subhi, tagint *tag,
                           int **nspecial, tagint **special,
                           int *nspecial15, tagint **special15,
                           const bool eflag_in, const bool vflag_in,
                           const bool eatom, const bool vatom, int &host_start,
                           int **ilist, int **jnum, const double cpu_time,
                           bool &success, double *host_q, double *boxlo,
                           double *prd, void **tep_ptr) {
  acc_timers();
  int eflag, vflag;
  if (eatom) eflag=2;
  else if (eflag_in) eflag=1;
  else eflag=0;
  if (vatom) vflag=2;
  else if (vflag_in) vflag=1;
  else vflag=0;

  #ifdef LAL_NO_BLOCK_REDUCE
  if (eflag) eflag=2;
  if (vflag) vflag=2;
  #endif

  set_kernel(eflag,vflag);

  // ------------------- Resize _tep array ------------------------

  if (nall>_max_tep_size) {
    _max_tep_size=static_cast<int>(static_cast<double>(nall)*1.10);
    _tep.resize(_max_tep_size*4);

    dev_nspecial15.clear();
    dev_special15.clear();
    dev_special15_t.clear();
    dev_nspecial15.alloc(nall,*(this->ucl_device),UCL_READ_ONLY);
    dev_special15.alloc(_maxspecial15*nall,*(this->ucl_device),UCL_READ_ONLY);
    dev_special15_t.alloc(nall*_maxspecial15,*(this->ucl_device),UCL_READ_ONLY);   
  }
  *tep_ptr=_tep.host.begin();

  if (inum_full==0) {
    host_start=0;
    // Make sure textures are correct if realloc by a different hybrid style
    resize_atom(0,nall,success);
    zero_timers();
    return nullptr;
  }

  hd_balancer.balance(cpu_time);
  int inum=hd_balancer.get_gpu_count(ago,inum_full);
  ans->inum(inum);
  host_start=inum;

  // Build neighbor list on GPU if necessary
  if (ago==0) {
    build_nbor_list(inum, inum_full-inum, nall, host_x, host_type,
                    sublo, subhi, tag, nspecial, special, nspecial15, special15,
                    success);
    if (!success)
      return nullptr;
    atom->cast_q_data(host_q);
    cast_extra_data(host_amtype, host_amgroup, host_rpole, host_uind, host_uinp);
    hd_balancer.start_timer();
  } else {
    atom->cast_x_data(host_x,host_type);
    atom->cast_q_data(host_q);
    cast_extra_data(host_amtype, host_amgroup, host_rpole, host_uind, host_uinp);
    hd_balancer.start_timer();
    atom->add_x_data(host_x,host_type);
  }
  atom->add_q_data();
  atom->add_extra_data();

  *ilist=nbor->host_ilist.begin();
  *jnum=nbor->host_acc.begin();

  device->precompute(ago,inum_full,nall,host_x,host_type,success,host_q,
                     boxlo, prd);

  const int red_blocks=loop(eflag,vflag);
  ans->copy_answers(eflag_in,vflag_in,eatom,vatom,red_blocks);
  device->add_ans_object(ans);
  hd_balancer.stop_timer();

  // copy tep from device to host

  _tep.update_host(_max_tep_size*4,false);
/*  
  printf("GPU lib: tep size = %d: max tep size = %d\n", this->_tep.cols(), _max_tep_size);
  for (int i = 0; i < 10; i++) {
    numtyp4* p = (numtyp4*)(&this->_tep[4*i]);
    printf("i = %d; tep = %f %f %f\n", i, p->x, p->y, p->z);
  }
*/  
  return nbor->host_jlist.begin()-host_start;
}

template <class numtyp, class acctyp>
double BaseAmoebaT::host_memory_usage_atomic() const {
  return device->atom.host_memory_usage()+nbor->host_memory_usage()+
         4*sizeof(numtyp)+sizeof(BaseAmoeba<numtyp,acctyp>);
}

template <class numtyp, class acctyp>
void BaseAmoebaT::cast_extra_data(int* amtype, int* amgroup, double** rpole,
    double** uind, double** uinp) {
  int _nall=atom->nall();
  numtyp *pextra=reinterpret_cast<numtyp*>(&(atom->extra[0]));

  int n = 0;
  int nstride = 4;
  for (int i = 0; i < _nall; i++) {
    int idx = n+i*nstride;
    pextra[idx]   = rpole[i][0];
    pextra[idx+1] = rpole[i][1];
    pextra[idx+2] = rpole[i][2];
    pextra[idx+3] = rpole[i][3];
  }

  n += nstride*_nall;
  for (int i = 0; i < _nall; i++) {
    int idx = n+i*nstride;
    pextra[idx]   = rpole[i][4];
    pextra[idx+1] = rpole[i][5];
    pextra[idx+2] = rpole[i][6];
    pextra[idx+3] = rpole[i][8];
  }

  n += nstride*_nall;
  for (int i = 0; i < _nall; i++) {
    int idx = n+i*nstride;
    pextra[idx]   = rpole[i][9];
    pextra[idx+1] = rpole[i][12];
    pextra[idx+2] = (numtyp)amtype[i];
    pextra[idx+3] = (numtyp)amgroup[i];
  }

  n += nstride*_nall;
  for (int i = 0; i < _nall; i++) {
    int idx = n+i*nstride;
    pextra[idx]   = uind[i][0];
    pextra[idx+1] = uind[i][1];
    pextra[idx+2] = uind[i][2];
  }

  n += nstride*_nall;
  for (int i = 0; i < _nall; i++) {
    int idx = n+i*nstride;
    pextra[idx]   = uinp[i][0];
    pextra[idx+1] = uinp[i][1];
    pextra[idx+2] = uinp[i][2];    
  }
}

template <class numtyp, class acctyp>
void BaseAmoebaT::compile_kernels(UCL_Device &dev, const void *pair_str,
                                  const char *kname) {
  if (_compiled)
    return;

  if (pair_program) delete pair_program;
  pair_program=new UCL_Program(dev);
  std::string oclstring = device->compile_string()+" -DEVFLAG=1";
  pair_program->load_string(pair_str,oclstring.c_str(),nullptr,screen);
  
  k_polar.set_function(*pair_program,kname);
  k_special15.set_function(*pair_program,"k_special15");
  pos_tex.get_texture(*pair_program,"pos_tex");
  q_tex.get_texture(*pair_program,"q_tex");
  
  _compiled=true;

  #if defined(USE_OPENCL) && (defined(CL_VERSION_2_1) || defined(CL_VERSION_3_0))
  if (dev.has_subgroup_support()) {
    size_t mx_subgroup_sz = k_polar.max_subgroup_size(_block_size);
    if (_threads_per_atom > mx_subgroup_sz)
      _threads_per_atom = mx_subgroup_sz;
    device->set_simd_size(mx_subgroup_sz);
  }
  #endif

}

template <class numtyp, class acctyp>
int BaseAmoebaT::add_onefive_neighbors() {
  // Compute the block size and grid size to keep all cores busy
  const int BX=block_size();
  int GX=static_cast<int>(ceil(static_cast<double>(ans->inum())/
                               (BX/_threads_per_atom)));

  int _nall=atom->nall();
  int ainum=ans->inum();
  int nbor_pitch=nbor->nbor_pitch();
  
  k_special15.set_size(GX,BX);
  k_special15.run(&nbor->dev_nbor, &_nbor_data->begin(),
                    &atom->dev_tag, &dev_nspecial15, &dev_special15,
                    &ainum, &_nall, &nbor_pitch,
                    &_threads_per_atom);
  
  return GX;
}

template class BaseAmoeba<PRECISION,ACC_PRECISION>;
}
