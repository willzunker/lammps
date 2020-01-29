/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#define MAX(A,B) ((A) > (B) ? (A) : (B))

// data types for 2d/3d FFTs

#ifndef FFT_DATA_KOKKOS_H
#define FFT_DATA_KOKKOS_H

// User-settable FFT precision

// FFT_PRECISION = 1 is single-precision complex (4-byte real, 4-byte imag)
// FFT_PRECISION = 2 is double-precision complex (8-byte real, 8-byte imag)

#ifdef FFT_SINGLE
#define FFT_PRECISION 1
#define MPI_FFT_SCALAR MPI_FLOAT
typedef float FFT_SCALAR;
#else
#define FFT_PRECISION 2
#define MPI_FFT_SCALAR MPI_DOUBLE
typedef double FFT_SCALAR;
#endif

// -------------------------------------------------------------------------

// Data types for single-precision complex

#if FFT_PRECISION == 1
#elif FFT_PRECISION == 2
#else
#error "FFT_PRECISION needs to be either 1 (=single) or 2 (=double)"
#endif


// with KOKKOS in CUDA mode we can only have
// CUFFT or KISSFFT, thus undefine all other
// FFTs here, since they may be valid in fft3d.cpp

#if defined(KOKKOS_ENABLE_CUDA)
# if defined(FFT_FFTW)
#  undef FFT_FFTW
# endif
# if defined(FFT_FFTW3)
#  undef FFT_FFTW3
# endif
# if defined(FFT_MKL)
#  undef FFT_MKL
# endif
# if !defined(FFT_CUFFT) && !defined(FFT_KISSFFT)
#  define FFT_KISSFFT
# endif
#else
# if defined(FFT_CUFFT)
#  error "Must enable CUDA with KOKKOS to use -DFFT_CUFFT"
# endif
// if user set FFTW, it means FFTW3
# ifdef FFT_FFTW
#  define FFT_FFTW3
# endif
# ifdef FFT_FFTW_THREADS
#  if !defined(FFT_FFTW3)
#   error "Must use -DFFT_FFTW3 with -DFFT_FFTW_THREADS"
#  endif
# endif
#endif

#if defined(FFT_MKL)
  #include "mkl_dfti.h"
  #if defined(FFT_SINGLE)
    typedef float _Complex FFT_DATA;
    #define FFT_MKL_PREC DFTI_SINGLE
  #else
    typedef double _Complex FFT_DATA;
    #define FFT_MKL_PREC DFTI_DOUBLE
  #endif
#elif defined(FFT_FFTW3)
  #include "fftw3.h"
  #if defined(FFT_SINGLE)
    typedef fftwf_complex FFT_DATA;
    #define FFTW_API(function)  fftwf_ ## function
  #else
    typedef fftw_complex FFT_DATA;
    #define FFTW_API(function) fftw_ ## function
  #endif
#elif defined(FFT_CUFFT)
  #include "cufft.h"
  #if defined(FFT_SINGLE)
    #define cufftExec cufftExecC2C
    #define CUFFT_TYPE CUFFT_C2C
    typedef cufftComplex FFT_DATA;
  #else
    #define cufftExec cufftExecZ2Z
    #define CUFFT_TYPE CUFFT_Z2Z
    typedef cufftDoubleComplex FFT_DATA;
  #endif
#else
  #include "kissfft_kokkos.h"
  #if defined(FFT_SINGLE)
    #define kiss_fft_scalar float
  #else
    #define kiss_fft_scalar double
    typedef struct {
        kiss_fft_scalar re;
        kiss_fft_scalar im;
    } FFT_DATA;
  #endif
  #ifndef FFT_KISSFFT
  #define FFT_KISSFFT
  #endif
#endif

#endif
