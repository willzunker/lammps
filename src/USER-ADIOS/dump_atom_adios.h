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

#ifdef DUMP_CLASS

DumpStyle(atom/adios,DumpAtomADIOS)

#else

#ifndef LMP_DUMP_ATOM_ADIOS_H
#define LMP_DUMP_ATOM_ADIOS_H

#include "dump_atom.h"
#include <stdlib.h>
#include <stdint.h>
#include "adios2.h"

namespace LAMMPS_NS {

class DumpAtomADIOS : public DumpAtom {

 public:
  DumpAtomADIOS(class LAMMPS *, int, char **);
  virtual ~DumpAtomADIOS();

 protected:

  const std::string ioName="atom";   // name of adios group, referrable in adios2_config.xml
  adios2::ADIOS *ad = nullptr; // adios object
  adios2::IO io;    // adios group of variables and attributes in this dump
  adios2::Engine fh; // adios file/stream handle object
  adios2::Variable<double> varAtoms; // one ADIOS output variable we need to change
  uint64_t groupSize; // pre-calculate # of bytes written per processor in a step before writing anything
  uint64_t groupTotalSize; // ADIOS buffer size returned by adios_group_size(), valid only if size is > default 16MB ADIOS buffer
  std::string filecurrent;  // name of file for this round (with % and * replaced)

  virtual void openfile();
  virtual void write();
  virtual void init_style();

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Cannot open dump file %s

The output file for the dump command cannot be opened.  Check that the
path and name are correct.

E: Too much per-proc info for dump

Number of local atoms times number of columns must fit in a 32-bit
integer for dump.

*/
