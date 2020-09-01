# ----------------------------------------------------------------------
#   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
#   http://lammps.sandia.gov, Sandia National Laboratories
#   Steve Plimpton, sjplimp@sandia.gov
#
#   Copyright (2003) Sandia Corporation.  Under the terms of Contract
#   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
#   certain rights in this software.  This software is distributed under
#   the GNU General Public License.
#
#   See the README file in the top-level LAMMPS directory.
# -------------------------------------------------------------------------
# Python wrappers for the LAMMPS library via ctypes

# for python2/3 compatibility

from __future__ import print_function

# imports for simple LAMMPS python wrapper module "lammps"

import sys,traceback,types
from ctypes import *
from os.path import dirname,abspath,join
from inspect import getsourcefile

# imports for advanced LAMMPS python wrapper modules "PyLammps" and "IPyLammps"

from collections import namedtuple
import os
import select
import re
import sys

# various symbolic constants to be used
# in certain calls to select data formats
LAMMPS_INT    = 0
LAMMPS_INT2D  = 1
LAMMPS_DOUBLE = 2
LAMMPS_DOUBLE2D = 3
LAMMPS_BIGINT = 4
LAMMPS_TAGINT = 5
LAMMPS_STRING = 6

# these must be kept in sync with the enums in library.h
LMP_STYLE_GLOBAL = 0
LMP_STYLE_ATOM   = 1
LMP_STYLE_LOCAL  = 2

LMP_TYPE_SCALAR  = 0
LMP_TYPE_VECTOR  = 1
LMP_TYPE_ARRAY   = 2
LMP_SIZE_VECTOR  = 3
LMP_SIZE_ROWS    = 4
LMP_SIZE_COLS    = 5

LMP_VAR_EQUAL = 0
LMP_VAR_ATOM  = 1

# -------------------------------------------------------------------------

def get_ctypes_int(size):
  if size == 4:
    return c_int32
  elif size == 8:
    return c_int64
  return c_int

# -------------------------------------------------------------------------

class MPIAbortException(Exception):
  def __init__(self, message):
    self.message = message

  def __str__(self):
    return repr(self.message)

# -------------------------------------------------------------------------

class NeighList:
    """This is a wrapper class that exposes the contents of a neighbor list.

    It can be used like a regular Python list.

    Internally it uses the lower-level LAMMPS C-library interface.

    :param lmp: reference to instance of :py:class:`lammps`
    :type  lmp: lammps
    :param idx: neighbor list index
    :type  idx: int
    """
    def __init__(self, lmp, idx):
        self.lmp = lmp
        self.idx = idx

    def __str__(self):
        return "Neighbor List ({} atoms)".format(self.size)

    def __repr__(self):
        return self.__str__()

    @property
    def size(self):
        """
        :return: number of elements in neighbor list
        """
        return self.lmp.get_neighlist_size(self.idx)

    def get(self, element):
        """
        :return: tuple with atom local index, number of neighbors and array of neighbor local atom indices
        :rtype:  (int, int, numpy.array)
        """
        iatom, numneigh, neighbors = self.lmp.get_neighlist_element_neighbors(self.idx, element)
        return iatom, numneigh, neighbors

    # the methods below implement the iterator interface, so NeighList can be used like a regular Python list

    def __getitem__(self, element):
        return self.get(element)

    def __len__(self):
        return self.size

    def __iter__(self):
        inum = self.size

        for ii in range(inum):
            yield self.get(ii)

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

class lammps(object):
  """Create an instance of the LAMMPS Python class.

  .. _mpi4py_docs: https://mpi4py.readthedocs.io/

  This is a Python wrapper class that exposes the LAMMPS C-library
  interface to Python.  It either requires that LAMMPS has been compiled
  as shared library which is then dynamically loaded via the ctypes
  Python module or that this module called from a Python function that
  is called from a Python interpreter embedded into a LAMMPS executable,
  for example through the :doc:`python invoke <python>` command.
  When the class is instantiated it calls the :cpp:func:`lammps_open`
  function of the LAMMPS C-library interface, which in
  turn will create an instance of the :cpp:class:`LAMMPS <LAMMPS_NS::LAMMPS>`
  C++ class.  The handle to this C++ class is stored internally
  and automatically passed to the calls to the C library interface.

  :param name: "machine" name of the shared LAMMPS library ("mpi" loads ``liblammps_mpi.so``, "" loads ``liblammps.so``)
  :type  name: string
  :param cmdargs: list of command line arguments to be passed to the :cpp:func:`lammps_open` function.  The executable name is automatically added.
  :type  cmdargs: list
  :param ptr: pointer to a LAMMPS C++ class instance when called from an embedded Python interpreter.  None means load symbols from shared library.
  :type  ptr: pointer
  :param comm: MPI communicator (as provided by `mpi4py <mpi4py_docs_>`_). ``None`` means use ``MPI_COMM_WORLD`` implicitly.
  :type  comm: MPI_Comm
  """

  # -------------------------------------------------------------------------
  # create an instance of LAMMPS

  def __init__(self,name='',cmdargs=None,ptr=None,comm=None):
    self.comm = comm
    self.opened = 0

    # determine module file location

    modpath = dirname(abspath(getsourcefile(lambda:0)))
    self.lib = None
    self.lmp = None

    # if a pointer to a LAMMPS object is handed in
    # when being called from a Python interpreter
    # embedded into a LAMMPS executable, all library
    # symbols should already be available so we do not
    # load a shared object.

    try:
      if ptr: self.lib = CDLL("",RTLD_GLOBAL)
    except:
      self.lib = None

    # load liblammps.so unless name is given
    #   if name = "g++", load liblammps_g++.so
    # try loading the LAMMPS shared object from the location
    #   of lammps.py with an absolute path,
    #   so that LD_LIBRARY_PATH does not need to be set for regular install
    # fall back to loading with a relative path,
    #   typically requires LD_LIBRARY_PATH to be set appropriately

    if any([f.startswith('liblammps') and f.endswith('.dylib')
            for f in os.listdir(modpath)]):
      lib_ext = ".dylib"
    elif any([f.startswith('liblammps') and f.endswith('.dll')
              for f in os.listdir(modpath)]):
      lib_ext = ".dll"
    else:
      lib_ext = ".so"

    if not self.lib:
      try:
        if not name:
          self.lib = CDLL(join(modpath,"liblammps" + lib_ext),RTLD_GLOBAL)
        else:
          self.lib = CDLL(join(modpath,"liblammps_%s" % name + lib_ext),
                          RTLD_GLOBAL)
      except:
        if not name:
          self.lib = CDLL("liblammps" + lib_ext,RTLD_GLOBAL)
        else:
          self.lib = CDLL("liblammps_%s" % name + lib_ext,RTLD_GLOBAL)


    # declare all argument and return types for all library methods here.
    # exceptions are where the arguments depend on certain conditions and
    # then are defined where the functions are used.
    self.lib.lammps_extract_setting.argtypes = [c_void_p, c_char_p]
    self.lib.lammps_extract_setting.restype = c_int

    # set default types
    # needed in later declarations
    self.c_bigint = get_ctypes_int(self.extract_setting("bigint"))
    self.c_tagint = get_ctypes_int(self.extract_setting("tagint"))
    self.c_imageint = get_ctypes_int(self.extract_setting("imageint"))

    self.lib.lammps_open.restype = c_void_p
    self.lib.lammps_open_no_mpi.restype = c_void_p
    self.lib.lammps_close.argtypes = [c_void_p]
    self.lib.lammps_free.argtypes = [c_void_p]

    self.lib.lammps_file.argtypes = [c_void_p, c_char_p]
    self.lib.lammps_file.restype = None

    self.lib.lammps_command.argtypes = [c_void_p, c_char_p]
    self.lib.lammps_command.restype = c_char_p
    self.lib.lammps_commands_list.restype = None
    self.lib.lammps_commands_string.argtypes = [c_void_p, c_char_p]
    self.lib.lammps_commands_string.restype = None

    self.lib.lammps_get_natoms.argtypes = [c_void_p]
    self.lib.lammps_get_natoms.restype = c_double
    self.lib.lammps_extract_box.argtypes = \
      [c_void_p,POINTER(c_double),POINTER(c_double),
       POINTER(c_double),POINTER(c_double),POINTER(c_double),
       POINTER(c_int),POINTER(c_int)]
    self.lib.lammps_extract_box.restype = None

    self.lib.lammps_reset_box.argtypes = \
      [c_void_p,POINTER(c_double),POINTER(c_double),c_double,c_double,c_double]
    self.lib.lammps_reset_box.restype = None

    self.lib.lammps_gather_atoms.argtypes = \
      [c_void_p,c_char_p,c_int,c_int,c_void_p]
    self.lib.lammps_gather_atoms.restype = None

    self.lib.lammps_gather_atoms_concat.argtypes = \
      [c_void_p,c_char_p,c_int,c_int,c_void_p]
    self.lib.lammps_gather_atoms_concat.restype = None

    self.lib.lammps_gather_atoms_subset.argtypes = \
      [c_void_p,c_char_p,c_int,c_int,c_int,POINTER(c_int),c_void_p]
    self.lib.lammps_gather_atoms_subset.restype = None

    self.lib.lammps_scatter_atoms.argtypes = \
      [c_void_p,c_char_p,c_int,c_int,c_void_p]
    self.lib.lammps_scatter_atoms.restype = None

    self.lib.lammps_scatter_atoms_subset.argtypes = \
      [c_void_p,c_char_p,c_int,c_int,c_int,POINTER(c_int),c_void_p]
    self.lib.lammps_scatter_atoms_subset.restype = None

    self.lib.lammps_gather.argtypes = \
      [c_void_p,c_char_p,c_int,c_int,c_void_p]
    self.lib.lammps_gather.restype = None

    self.lib.lammps_gather_concat.argtypes = \
      [c_void_p,c_char_p,c_int,c_int,c_void_p]
    self.lib.lammps_gather_concat.restype = None

    self.lib.lammps_gather_subset.argtypes = \
      [c_void_p,c_char_p,c_int,c_int,c_int,POINTER(c_int),c_void_p]
    self.lib.lammps_gather_subset.restype = None

    self.lib.lammps_scatter.argtypes = \
      [c_void_p,c_char_p,c_int,c_int,c_void_p]
    self.lib.lammps_scatter.restype = None

    self.lib.lammps_scatter_subset.argtypes = \
      [c_void_p,c_char_p,c_int,c_int,c_int,POINTER(c_int),c_void_p]
    self.lib.lammps_scatter_subset.restype = None


    self.lib.lammps_find_pair_neighlist.argtypes = [c_void_p, c_char_p, c_int, c_int, c_int]
    self.lib.lammps_find_pair_neighlist.restype  = c_int

    self.lib.lammps_find_fix_neighlist.argtypes = [c_void_p, c_char_p, c_int]
    self.lib.lammps_find_fix_neighlist.restype  = c_int

    self.lib.lammps_find_compute_neighlist.argtypes = [c_void_p, c_char_p, c_int]
    self.lib.lammps_find_compute_neighlist.restype  = c_int

    self.lib.lammps_neighlist_num_elements.argtypes = [c_void_p, c_int]
    self.lib.lammps_neighlist_num_elements.restype  = c_int

    self.lib.lammps_neighlist_element_neighbors.argtypes = [c_void_p, c_int, c_int, POINTER(c_int), POINTER(c_int), POINTER(POINTER(c_int))]
    self.lib.lammps_neighlist_element_neighbors.restype  = None

    self.lib.lammps_has_error.argtypes = [c_void_p]
    self.lib.lammps_has_error.restype = c_bool

    self.lib.lammps_get_last_error_message.argtypes = [c_void_p, c_char_p, c_int]
    self.lib.lammps_get_last_error_message.restype = c_int

    self.lib.lammps_extract_global.argtypes = [c_void_p, c_char_p]
    self.lib.lammps_extract_compute.argtypes = [c_void_p, c_char_p, c_int, c_int]

    self.lib.lammps_get_thermo.argtypes = [c_void_p, c_char_p]
    self.lib.lammps_get_thermo.restype = c_double

    self.lib.lammps_encode_image_flags.restype = self.c_imageint

    self.lib.lammps_config_package_name.argtypes = [c_int, c_char_p, c_int]

    self.lib.lammps_has_style.argtypes = [c_void_p, c_char_p, c_char_p]

    self.lib.lammps_set_variable.argtypes = [c_void_p, c_char_p, c_char_p]

    self.lib.lammps_style_count.argtypes = [c_void_p, c_char_p]

    self.lib.lammps_style_name.argtypes = [c_void_p, c_char_p, c_int, c_char_p, c_int]

    self.lib.lammps_version.argtypes = [c_void_p]

    self.lib.lammps_decode_image_flags.argtypes = [self.c_imageint, POINTER(c_int*3)]

    self.lib.lammps_extract_atom.argtypes = [c_void_p, c_char_p]

    self.lib.lammps_extract_fix.argtypes = [c_void_p, c_char_p, c_int, c_int, c_int, c_int]

    self.lib.lammps_extract_variable.argtypes = [c_void_p, c_char_p, c_char_p]

    # TODO: NOT IMPLEMENTED IN PYTHON WRAPPER
    self.lib.lammps_fix_external_set_energy_global = [c_void_p, c_char_p, c_double]
    self.lib.lammps_fix_external_set_virial_global = [c_void_p, c_char_p, POINTER(c_double)]

    # detect if Python is using version of mpi4py that can pass a communicator

    self.has_mpi4py = False
    try:
      from mpi4py import __version__ as mpi4py_version
      # tested to work with mpi4py versions 2 and 3
      self.has_mpi4py = mpi4py_version.split('.')[0] in ['2','3']
    except:
      pass

    # if no ptr provided, create an instance of LAMMPS
    #   don't know how to pass an MPI communicator from PyPar
    #   but we can pass an MPI communicator from mpi4py v2.0.0 and later
    #   no_mpi call lets LAMMPS use MPI_COMM_WORLD
    #   cargs = array of C strings from args
    # if ptr, then are embedding Python in LAMMPS input script
    #   ptr is the desired instance of LAMMPS
    #   just convert it to ctypes ptr and store in self.lmp

    if not ptr:

      # with mpi4py v2, can pass MPI communicator to LAMMPS
      # need to adjust for type of MPI communicator object
      # allow for int (like MPICH) or void* (like OpenMPI)
      if self.has_mpi4py and self.has_mpi_support:
        from mpi4py import MPI
        self.MPI = MPI

      if comm:
        if not self.has_mpi4py:
          raise Exception('Python mpi4py version is not 2 or 3')
        if not self.has_mpi_support:
          raise Exception('LAMMPS not compiled with real MPI library')
        if self.MPI._sizeof(self.MPI.Comm) == sizeof(c_int):
          MPI_Comm = c_int
        else:
          MPI_Comm = c_void_p

        narg = 0
        cargs = None
        if cmdargs:
          cmdargs.insert(0,"lammps.py")
          narg = len(cmdargs)
          for i in range(narg):
            if type(cmdargs[i]) is str:
              cmdargs[i] = cmdargs[i].encode()
          cargs = (c_char_p*narg)(*cmdargs)
          self.lib.lammps_open.argtypes = [c_int, c_char_p*narg, \
                                           MPI_Comm, c_void_p]
        else:
          self.lib.lammps_open.argtypes = [c_int, c_char_p, \
                                           MPI_Comm, c_void_p]

        self.opened = 1
        comm_ptr = self.MPI._addressof(comm)
        comm_val = MPI_Comm.from_address(comm_ptr)
        self.lmp = c_void_p(self.lib.lammps_open(narg,cargs,comm_val,None))

      else:
        if self.has_mpi4py and self.has_mpi_support:
          self.comm = self.MPI.COMM_WORLD
        self.opened = 1
        if cmdargs:
          cmdargs.insert(0,"lammps.py")
          narg = len(cmdargs)
          for i in range(narg):
            if type(cmdargs[i]) is str:
              cmdargs[i] = cmdargs[i].encode()
          cargs = (c_char_p*narg)(*cmdargs)
          self.lib.lammps_open_no_mpi.argtypes = [c_int, c_char_p*narg, \
                                                  c_void_p]
          self.lmp = c_void_p(self.lib.lammps_open_no_mpi(narg,cargs,None))
        else:
          self.lib.lammps_open_no_mpi.argtypes = [c_int, c_char_p, c_void_p]
          self.lmp = c_void_p(self.lib.lammps_open_no_mpi(0,None,None))

    else:
      # magic to convert ptr to ctypes ptr
      if sys.version_info >= (3, 0):
        # Python 3 (uses PyCapsule API)
        pythonapi.PyCapsule_GetPointer.restype = c_void_p
        pythonapi.PyCapsule_GetPointer.argtypes = [py_object, c_char_p]
        self.lmp = c_void_p(pythonapi.PyCapsule_GetPointer(ptr, None))
      else:
        # Python 2 (uses PyCObject API)
        pythonapi.PyCObject_AsVoidPtr.restype = c_void_p
        pythonapi.PyCObject_AsVoidPtr.argtypes = [py_object]
        self.lmp = c_void_p(pythonapi.PyCObject_AsVoidPtr(ptr))

    # optional numpy support (lazy loading)
    self._numpy = None

    self._installed_packages = None
    self._available_styles = None

    # add way to insert Python callback for fix external
    self.callback = {}
    self.FIX_EXTERNAL_CALLBACK_FUNC = CFUNCTYPE(None, py_object, self.c_bigint, c_int, POINTER(self.c_tagint), POINTER(POINTER(c_double)), POINTER(POINTER(c_double)))
    self.lib.lammps_set_fix_external_callback.argtypes = [c_void_p, c_char_p, self.FIX_EXTERNAL_CALLBACK_FUNC, py_object]
    self.lib.lammps_set_fix_external_callback.restype = None

  # -------------------------------------------------------------------------
  # shut-down LAMMPS instance

  def __del__(self):
    if self.lmp and self.opened:
      self.lib.lammps_close(self.lmp)
      self.opened = 0

  # -------------------------------------------------------------------------

  @property
  def numpy(self):
    "Convert between ctypes arrays and numpy arrays"
    if not self._numpy:
      import numpy as np
      class LammpsNumpyWrapper:
        def __init__(self, lmp):
          self.lmp = lmp

        def _ctype_to_numpy_int(self, ctype_int):
          if ctype_int == c_int32:
            return np.int32
          elif ctype_int == c_int64:
            return np.int64
          return np.intc

        def extract_atom_iarray(self, name, nelem, dim=1):
          if name in ['id', 'molecule']:
            c_int_type = self.lmp.c_tagint
          elif name in ['image']:
            c_int_type = self.lmp.c_imageint
          else:
            c_int_type = c_int

          if dim == 1:
            raw_ptr = self.lmp.extract_atom(name, LAMMPS_INT)
          else:
            raw_ptr = self.lmp.extract_atom(name, LAMMPS_INT2D)

          return self.iarray(c_int_type, raw_ptr, nelem, dim)

        def extract_atom_darray(self, name, nelem, dim=1):
          if dim == 1:
            raw_ptr = self.lmp.extract_atom(name, LAMMPS_DOUBLE)
          else:
            raw_ptr = self.lmp.extract_atom(name, LAMMPS_DOUBLE2D)

          return self.darray(raw_ptr, nelem, dim)

        def extract_compute(self, cid, style, datatype):
          value = self.lmp.extract_compute(cid, style, datatype)

          if style in (LMP_STYLE_GLOBAL, LMP_STYLE_LOCAL):
            if datatype == LMP_TYPE_VECTOR:
              nrows = self.lmp.extract_compute(cid, style, LMP_SIZE_VECTOR)
              return self.darray(value, nrows)
            elif datatype == LMP_TYPE_ARRAY:
              nrows = self.lmp.extract_compute(cid, style, LMP_SIZE_ROWS)
              ncols = self.lmp.extract_compute(cid, style, LMP_SIZE_COLS)
              return self.darray(value, nrows, ncols)
          elif style == LMP_STYLE_ATOM:
            if datatype == LMP_TYPE_VECTOR:
              nlocal = self.lmp.extract_global("nlocal", LAMMPS_INT)
              return self.darray(value, nlocal)
            elif datatype == LMP_TYPE_ARRAY:
              nlocal = self.lmp.extract_global("nlocal", LAMMPS_INT)
              ncols = self.lmp.extract_compute(cid, style, LMP_SIZE_COLS)
              return self.darray(value, nlocal, ncols)
          return value

        def extract_fix(self, fid, style, datatype, nrow=0, ncol=0):
          value = self.lmp.extract_fix(fid, style, datatype, nrow, ncol)
          if style == LMP_STYLE_ATOM:
            if datatype == LMP_TYPE_VECTOR:
              nlocal = self.lmp.extract_global("nlocal", LAMMPS_INT)
              return self.darray(value, nlocal)
            elif datatype == LMP_TYPE_ARRAY:
              nlocal = self.lmp.extract_global("nlocal", LAMMPS_INT)
              ncols = self.lmp.extract_fix(fid, style, LMP_SIZE_COLS, 0, 0)
              return self.darray(value, nlocal, ncols)
          elif style == LMP_STYLE_LOCAL:
            if datatype == LMP_TYPE_VECTOR:
              nrows = self.lmp.extract_fix(fid, style, LMP_SIZE_ROWS, 0, 0)
              return self.darray(value, nrows)
            elif datatype == LMP_TYPE_ARRAY:
              nrows = self.lmp.extract_fix(fid, style, LMP_SIZE_ROWS, 0, 0)
              ncols = self.lmp.extract_fix(fid, style, LMP_SIZE_COLS, 0, 0)
              return self.darray(value, nrows, ncols)
          return value

        def extract_variable(self, name, group=None, datatype=LMP_VAR_EQUAL):
          value = self.lmp.extract_variable(name, group, datatype)
          if datatype == LMP_VAR_ATOM:
            return np.ctypeslib.as_array(value)
          return value

        def iarray(self, c_int_type, raw_ptr, nelem, dim=1):
          np_int_type = self._ctype_to_numpy_int(c_int_type)

          if dim == 1:
            ptr = cast(raw_ptr, POINTER(c_int_type * nelem))
          else:
            ptr = cast(raw_ptr[0], POINTER(c_int_type * nelem * dim))

          a = np.frombuffer(ptr.contents, dtype=np_int_type)
          a.shape = (nelem, dim)
          return a

        def darray(self, raw_ptr, nelem, dim=1):
          if dim == 1:
            ptr = cast(raw_ptr, POINTER(c_double * nelem))
          else:
            ptr = cast(raw_ptr[0], POINTER(c_double * nelem * dim))

          a = np.frombuffer(ptr.contents)
          a.shape = (nelem, dim)
          return a

      self._numpy = LammpsNumpyWrapper(self)
    return self._numpy

  # -------------------------------------------------------------------------

  def close(self):
    """Explicitly delete a LAMMPS instance through the C-library interface.

    This is a wrapper around the :cpp:func:`lammps_close` function of the C-library interface.
    """
    if self.opened: self.lib.lammps_close(self.lmp)
    self.lmp = None
    self.opened = 0

  # -------------------------------------------------------------------------

  def finalize(self):
    """Shut down the MPI communication through the library interface by calling :cpp:func:`lammps_finalize`.
    """
    if self.opened: self.lib.lammps_close(self.lmp)
    self.lmp = None
    self.opened = 0
    self.lib.lammps_finalize()

  # -------------------------------------------------------------------------

  def version(self):
    """Return a numerical representation of the LAMMPS version in use.

    This is a wrapper around the :cpp:func:`lammps_close` function of the C-library interface.

    :return: version number
    :rtype:  int
    """
    return self.lib.lammps_version(self.lmp)

  # -------------------------------------------------------------------------

  def file(self, path):
    """Read LAMMPS commands from a file.

    This is a wrapper around the :cpp:func:`lammps_file` function of the C-library interface.
    It will open the file with the name/path `file` and process the LAMMPS commands line by line until
    the end. The function will return when the end of the file is reached.

    :param path: Name of the file/path with LAMMPS commands
    :type path:  string
    """
    if path: path = path.encode()
    else: return
    self.lib.lammps_file(self.lmp, path)

  # -------------------------------------------------------------------------

  def command(self,cmd):
    """Process a single LAMMPS input command from a string.

    This is a wrapper around the :cpp:func:`lammps_command`
    function of the C-library interface.

    :param cmd: a single lammps command
    :type cmd:  string
    """
    if cmd: cmd = cmd.encode()
    else: return
    self.lib.lammps_command(self.lmp,cmd)

    if self.has_exceptions and self.lib.lammps_has_error(self.lmp):
      sb = create_string_buffer(100)
      error_type = self.lib.lammps_get_last_error_message(self.lmp, sb, 100)
      error_msg = sb.value.decode().strip()

      if error_type == 2:
        raise MPIAbortException(error_msg)
      raise Exception(error_msg)

  # -------------------------------------------------------------------------

  def commands_list(self,cmdlist):
    """Process multiple LAMMPS input commands from a list of strings.

    This is a wrapper around the
    :cpp:func:`lammps_commands_list` function of
    the C-library interface.

    :param cmdlist: a single lammps command
    :type cmdlist:  list of strings
    """
    cmds = [x.encode() for x in cmdlist if type(x) is str]
    narg = len(cmdlist)
    args = (c_char_p * narg)(*cmds)
    self.lib.lammps_commands_list.argtypes = [c_void_p, c_int, c_char_p * narg]
    self.lib.lammps_commands_list(self.lmp,narg,args)

  # -------------------------------------------------------------------------

  def commands_string(self,multicmd):
    """Process a block of LAMMPS input commands from a string.

    This is a wrapper around the
    :cpp:func:`lammps_commands_string`
    function of the C-library interface.

    :param multicmd: text block of lammps commands
    :type multicmd:  string
    """
    if type(multicmd) is str: multicmd = multicmd.encode()
    self.lib.lammps_commands_string(self.lmp,c_char_p(multicmd))

  # -------------------------------------------------------------------------

  def get_natoms(self):
    """Get the total number of atoms in the LAMMPS instance.

    Will be precise up to 53-bit signed integer due to the
    underlying :cpp:func:`lammps_get_natoms` function returning a double.

    :return: number of atoms
    :rtype: float
    """
    return self.lib.lammps_get_natoms(self.lmp)

  # -------------------------------------------------------------------------

  def extract_box(self):
    """Extract simulation box parameters

    This is a wrapper around the :cpp:func:`lammps_extract_box` function
    of the C-library interface.  Unlike in the C function, the result is
    returned as a list.

    :return: list of the extracted data: boxlo, boxhi, xy, yz, xz, periodicity, box_change
    :rtype: [ 3*double, 3*double, double, double, 3*int, int]
    """
    boxlo = (3*c_double)()
    boxhi = (3*c_double)()
    xy = c_double()
    yz = c_double()
    xz = c_double()
    periodicity = (3*c_int)()
    box_change = c_int()

    self.lib.lammps_extract_box(self.lmp,boxlo,boxhi,
                                byref(xy),byref(yz),byref(xz),
                                periodicity,byref(box_change))

    boxlo = boxlo[:3]
    boxhi = boxhi[:3]
    xy = xy.value
    yz = yz.value
    xz = xz.value
    periodicity = periodicity[:3]
    box_change = box_change.value

    return boxlo,boxhi,xy,yz,xz,periodicity,box_change

  # -------------------------------------------------------------------------

  def reset_box(self,boxlo,boxhi,xy,yz,xz):
    """Reset simulation box parameters

    This is a wrapper around the :cpp:func:`lammps_reset_box` function
    of the C-library interface.

    :param boxlo: new lower box boundaries
    :type boxlo: list of 3 floating point numbers
    :param boxhi: new upper box boundaries
    :type boxhi: list of 3 floating point numbers
    :param xy: xy tilt factor
    :type xy: float
    :param yz: yz tilt factor
    :type yz: float
    :param xz: xz tilt factor
    :type xz: float
    """
    cboxlo = (3*c_double)(*boxlo)
    cboxhi = (3*c_double)(*boxhi)
    self.lib.lammps_reset_box(self.lmp,cboxlo,cboxhi,xy,yz,xz)

  # -------------------------------------------------------------------------

  def get_thermo(self,name):
    """Get current value of a thermo keyword

    This is a wrapper around the :cpp:func:`lammps_get_thermo`
    function of the C-library interface.

    :param name: name of thermo keyword
    :type name: string
    :return: value of thermo keyword
    :rtype: double or None
    """
    if name: name = name.encode()
    else: return None
    return self.lib.lammps_get_thermo(self.lmp,name)

  # -------------------------------------------------------------------------

  def extract_setting(self, name):
    """Query LAMMPS about global settings that can be expressed as an integer.

    This is a wrapper around the :cpp:func:`lammps_extract_setting`
    function of the C-library interface.  Its documentation includes
    a list of the supported keywords.

    :param name: name of the setting
    :type name:  string
    :return: value of the setting
    :rtype: int
    """
    if name: name = name.encode()
    else: return None
    return int(self.lib.lammps_extract_setting(self.lmp,name))

  # -------------------------------------------------------------------------
  # extract global info

  def extract_global(self, name, type):
    """Query LAMMPS about global settings of different types.

    This is a wrapper around the :cpp:func:`lammps_extract_global`
    function of the C-library interface.  Unlike the C function
    this method returns the value and not a pointer and thus can
    only return the first value for keywords representing a list
    of values.  The :cpp:func:`lammps_extract_global` documentation
    includes a list of the supported keywords and their data types.
    Since Python needs to know the data type to be able to interpret
    the result, the type has to be provided as an argument.  For
    that purpose the :py:mod:`lammps` module contains the constants
    ``LAMMPS_INT``, ``LAMMPS_DOUBLE``, ``LAMMPS_BIGINT``,
    ``LAMMPS_TAGINT``, and ``LAMMPS_STRING``.
    This function returns ``None`` if either the keyword is not
    recognized, or an invalid data type constant is used.

    :param name: name of the setting
    :type name:  string
    :param type: type of the returned data
    :type type:  int
    :return: value of the setting
    :rtype: integer or double or string or None
    """
    if name: name = name.encode()
    else: return None
    if type == LAMMPS_INT:
      self.lib.lammps_extract_global.restype = POINTER(c_int)
    elif type == LAMMPS_DOUBLE:
      self.lib.lammps_extract_global.restype = POINTER(c_double)
    elif type == LAMMPS_BIGINT:
      self.lib.lammps_extract_global.restype = POINTER(self.c_bigint)
    elif type == LAMMPS_TAGINT:
      self.lib.lammps_extract_global.restype = POINTER(self.c_tagint)
    elif type == LAMMPS_STRING:
      self.lib.lammps_extract_global.restype = c_char_p
      ptr = self.lib.lammps_extract_global(self.lmp,name)
      return str(ptr,'ascii')
    else: return None
    ptr = self.lib.lammps_extract_global(self.lmp,name)
    if ptr: return ptr[0]
    else: return None

  # -------------------------------------------------------------------------
  # extract per-atom info
  # NOTE: need to insure are converting to/from correct Python type
  #   e.g. for Python list or NumPy or ctypes

  def extract_atom(self,name,type):
    """Retrieve per-atom properties from LAMMPS

    This is a wrapper around the :cpp:func:`lammps_extract_atom`
    function of the C-library interface. Its documentation includes a
    list of the supported keywords and their data types.
    Since Python needs to know the data type to be able to interpret
    the result, the type has to be provided as an argument.  For
    that purpose the :py:mod:`lammps` module contains the constants
    ``LAMMPS_INT``, ``LAMMPS_INT2D``, ``LAMMPS_DOUBLE``,
    and ``LAMMPS_DOUBLE2D``.
    This function returns ``None`` if either the keyword is not
    recognized, or an invalid data type constant is used.

    .. note::

       While the returned arrays of per-atom data are dimensioned
       for the range [0:nmax] - as is the underlying storage -
       the data is usually only valid for the range of [0:nlocal],
       unless the property of interest is also updated for ghost
       atoms.  In some cases, this depends on a LAMMPS setting, see
       for example :doc:`comm_modify vel yes <comm_modify>`.

    :param name: name of the setting
    :type name:  string
    :param type: type of the returned data
    :type type:  int
    :return: requested data
    :rtype: pointer to integer or double or None
    """
    ntypes = int(self.extract_setting('ntypes'))
    nmax   = int(self.extract_setting('nmax'))
    if name: name = name.encode()
    else: return None
    if type == LAMMPS_INT:
      self.lib.lammps_extract_atom.restype = POINTER(c_int)
    elif type == LAMMPS_INT2D:
      self.lib.lammps_extract_atom.restype = POINTER(POINTER(c_int))
    elif type == LAMMPS_DOUBLE:
      self.lib.lammps_extract_atom.restype = POINTER(c_double)
    elif type == LAMMPS_DOUBLE2D:
      self.lib.lammps_extract_atom.restype = POINTER(POINTER(c_double))
    else: return None
    ptr = self.lib.lammps_extract_atom(self.lmp,name)
    if ptr: return ptr
    else:   return None


  # -------------------------------------------------------------------------

  def extract_compute(self,id,style,type):
    """Retrieve data from a LAMMPS compute

    This is a wrapper around the :cpp:func:`lammps_extract_compute`
    function of the C-library interface.
    This function returns ``None`` if either the compute id is not
    recognized, or an invalid combination of :ref:`style <py_style_constants>`
    and :ref:`type <py_type_constants>` constants is used. The
    names and functionality of the constants are the same as for
    the corresponding C-library function.  For requests to return
    a scalar or a size, the value is returned, otherwise a pointer.

    :param id: compute ID
    :type id:  string
    :param style: style of the data retrieve (global, atom, or local)
    :type style:  int
    :param type: type or size of the returned data (scalar, vector, or array)
    :type type:  int
    :return: requested data
    :rtype: integer or double or pointer to 1d or 2d double array or None
    """
    if id: id = id.encode()
    else: return None

    if type == LMP_TYPE_SCALAR:
      if style == LMP_STYLE_GLOBAL:
        self.lib.lammps_extract_compute.restype = POINTER(c_double)
        ptr = self.lib.lammps_extract_compute(self.lmp,id,style,type)
        return ptr[0]
      elif style == LMP_STYLE_ATOM:
        return None
      elif style == LMP_STYLE_LOCAL:
        self.lib.lammps_extract_compute.restype = POINTER(c_int)
        ptr = self.lib.lammps_extract_compute(self.lmp,id,style,type)
        return ptr[0]

    if type == LMP_TYPE_VECTOR:
      self.lib.lammps_extract_compute.restype = POINTER(c_double)
      ptr = self.lib.lammps_extract_compute(self.lmp,id,style,type)
      return ptr

    if type == LMP_TYPE_ARRAY:
      self.lib.lammps_extract_compute.restype = POINTER(POINTER(c_double))
      ptr = self.lib.lammps_extract_compute(self.lmp,id,style,type)
      return ptr

    if type == LMP_SIZE_COLS:
      if style == LMP_STYLE_GLOBAL  \
         or style == LMP_STYLE_ATOM \
         or style == LMP_STYLE_LOCAL:
        self.lib.lammps_extract_compute.restype = POINTER(c_int)
        ptr = self.lib.lammps_extract_compute(self.lmp,id,style,type)
        return ptr[0]

    if type == LMP_SIZE_VECTOR  \
       or type == LMP_SIZE_ROWS:
      if style == LMP_STYLE_GLOBAL  \
         or style == LMP_STYLE_LOCAL:
        self.lib.lammps_extract_compute.restype = POINTER(c_int)
        ptr = self.lib.lammps_extract_compute(self.lmp,id,style,type)
        return ptr[0]

    return None

  # -------------------------------------------------------------------------
  # extract fix info
  # in case of global data, free memory for 1 double via lammps_free()
  # double was allocated by library interface function

  def extract_fix(self,id,style,type,nrow=0,ncol=0):
    """Retrieve data from a LAMMPS fix

    This is a wrapper around the :cpp:func:`lammps_extract_fix`
    function of the C-library interface.
    This function returns ``None`` if either the fix id is not
    recognized, or an invalid combination of :ref:`style <py_style_constants>`
    and :ref:`type <py_type_constants>` constants is used. The
    names and functionality of the constants are the same as for
    the corresponding C-library function.  For requests to return
    a scalar or a size, the value is returned, also when accessing
    global vectors or arrays, otherwise a pointer.

    :param id: fix ID
    :type id:  string
    :param style: style of the data retrieve (global, atom, or local)
    :type style:  int
    :param type: type or size of the returned data (scalar, vector, or array)
    :type type:  int
    :param nrow: index of global vector element or row index of global array element
    :type nrow:  int
    :param ncol: column index of global array element
    :type ncol:  int
    :return: requested data
    :rtype: integer or double value, pointer to 1d or 2d double array  or None

    """
    if id: id = id.encode()
    else: return None

    if style == LMP_STYLE_GLOBAL:
      if type == LMP_TYPE_SCALAR    \
         or type == LMP_TYPE_VECTOR \
         or type == LMP_TYPE_ARRAY:
        self.lib.lammps_extract_fix.restype = POINTER(c_double)
        ptr = self.lib.lammps_extract_fix(self.lmp,id,style,type,nrow,ncol)
        result = ptr[0]
        self.lib.lammps_free(ptr)
        return result
      elif type == LMP_SIZE_VECTOR  \
           or type == LMP_SIZE_ROWS \
           or type == LMP_SIZE_COLS:
        self.lib.lammps_extract_fix.restype = POINTER(c_int)
        ptr = self.lib.lammps_extract_fix(self.lmp,id,style,type,nrow,ncol)
        return ptr[0]
      else:
        return None

    elif style == LMP_STYLE_ATOM:
      if type == LMP_TYPE_VECTOR:
        self.lib.lammps_extract_fix.restype = POINTER(c_double)
      elif type == LMP_TYPE_ARRAY:
        self.lib.lammps_extract_fix.restype = POINTER(POINTER(c_double))
      elif type == LMP_SIZE_COLS:
        self.lib.lammps_extract_fix.restype = POINTER(c_int)
      else:
        return None
      ptr = self.lib.lammps_extract_fix(self.lmp,id,style,type,nrow,ncol)
      if type == LMP_SIZE_COLS:
        return ptr[0]
      else:
        return ptr

    elif style == LMP_STYLE_LOCAL:
      if type == LMP_TYPE_VECTOR:
        self.lib.lammps_extract_fix.restype = POINTER(c_double)
      elif type == LMP_TYPE_ARRAY:
        self.lib.lammps_extract_fix.restype = POINTER(POINTER(c_double))
      elif type == LMP_TYPE_SCALAR    \
           or type == LMP_SIZE_VECTOR \
           or type == LMP_SIZE_ROWS   \
           or type == LMP_SIZE_COLS:
        self.lib.lammps_extract_fix.restype = POINTER(c_int)
      else:
        return None
      ptr = self.lib.lammps_extract_fix(self.lmp,id,style,type,nrow,ncol)
      if type == LMP_TYPE_VECTOR or type == LMP_TYPE_ARRAY:
        return ptr
      else:
        return ptr[0]
    else:
      return None

  # -------------------------------------------------------------------------
  # extract variable info
  # free memory for 1 double or 1 vector of doubles via lammps_free()
  # for vector, must copy nlocal returned values to local c_double vector
  # memory was allocated by library interface function

  def extract_variable(self,name,group=None,type=LMP_VAR_EQUAL):
    """ Evaluate a LAMMPS variable and return its data

    This function is a wrapper around the function
    :cpp:func:`lammps_extract_variable` of the C-library interface,
    evaluates variable name and returns a copy of the computed data.
    The memory temporarily allocated by the C-interface is deleted
    after the data is copied to a python variable or list.
    The variable must be either an equal-style (or equivalent)
    variable or an atom-style variable. The variable type has to
    provided as type parameter which may be two constants:
    ``LMP_VAR_EQUAL`` or ``LMP_VAR_STRING``; it defaults to
    equal-style variables.
    The group parameter is only used for atom-style variables and
    defaults to the group "all" if set to ``None``, which is the default.

    :param name: name of the variable to execute
    :type name: string
    :param group: name of group for atom style variable
    :type group: string
    :param type: type of variable
    :type type: int
    :return: the requested data
    :rtype: double, array of doubles, or None
    """
    if name: name = name.encode()
    else: return None
    if group: group = group.encode()
    if type == LMP_VAR_EQUAL:
      self.lib.lammps_extract_variable.restype = POINTER(c_double)
      ptr = self.lib.lammps_extract_variable(self.lmp,name,group)
      result = ptr[0]
      self.lib.lammps_free(ptr)
      return result
    if type == LMP_VAR_ATOM:
      self.lib.lammps_extract_global.restype = POINTER(c_int)
      nlocalptr = self.lib.lammps_extract_global(self.lmp,"nlocal".encode())
      nlocal = nlocalptr[0]
      result = (c_double*nlocal)()
      self.lib.lammps_extract_variable.restype = POINTER(c_double)
      ptr = self.lib.lammps_extract_variable(self.lmp,name,group)
      for i in range(nlocal): result[i] = ptr[i]
      self.lib.lammps_free(ptr)
      return result
    return None

  # -------------------------------------------------------------------------

  def set_variable(self,name,value):
    """Set a new value for a LAMMPS string style variable

    This is a wrapper around the :cpp:func:`lammps_set_variable`
    function of the C-library interface.

    :param name: name of the variable
    :type name: string
    :param value: new variable value
    :type value: any. will be converted to a string
    :return: either 0 on success or -1 on failure
    :rtype: int
    """
    if name: name = name.encode()
    else: return -1
    if value: value = str(value).encode()
    else: return -1
    return self.lib.lammps_set_variable(self.lmp,name,value)

  # -------------------------------------------------------------------------

  # return vector of atom properties gathered across procs
  # 3 variants to match src/library.cpp
  # name = atom property recognized by LAMMPS in atom->extract()
  # type = 0 for integer values, 1 for double values
  # count = number of per-atom valus, 1 for type or charge, 3 for x or f
  # returned data is a 1d vector - doc how it is ordered?
  # NOTE: need to insure are converting to/from correct Python type
  #   e.g. for Python list or NumPy or ctypes

  def gather_atoms(self,name,type,count):
    if name: name = name.encode()
    natoms = self.lib.lammps_get_natoms(self.lmp)
    if type == 0:
      data = ((count*natoms)*c_int)()
      self.lib.lammps_gather_atoms(self.lmp,name,type,count,data)
    elif type == 1:
      data = ((count*natoms)*c_double)()
      self.lib.lammps_gather_atoms(self.lmp,name,type,count,data)
    else: return None
    return data

  # -------------------------------------------------------------------------

  def gather_atoms_concat(self,name,type,count):
    if name: name = name.encode()
    natoms = self.lib.lammps_get_natoms(self.lmp)
    if type == 0:
      data = ((count*natoms)*c_int)()
      self.lib.lammps_gather_atoms_concat(self.lmp,name,type,count,data)
    elif type == 1:
      data = ((count*natoms)*c_double)()
      self.lib.lammps_gather_atoms_concat(self.lmp,name,type,count,data)
    else: return None
    return data

  def gather_atoms_subset(self,name,type,count,ndata,ids):
    if name: name = name.encode()
    if type == 0:
      data = ((count*ndata)*c_int)()
      self.lib.lammps_gather_atoms_subset(self.lmp,name,type,count,ndata,ids,data)
    elif type == 1:
      data = ((count*ndata)*c_double)()
      self.lib.lammps_gather_atoms_subset(self.lmp,name,type,count,ndata,ids,data)
    else: return None
    return data

  # -------------------------------------------------------------------------

  # scatter vector of atom properties across procs
  # 2 variants to match src/library.cpp
  # name = atom property recognized by LAMMPS in atom->extract()
  # type = 0 for integer values, 1 for double values
  # count = number of per-atom valus, 1 for type or charge, 3 for x or f
  # assume data is of correct type and length, as created by gather_atoms()
  # NOTE: need to insure are converting to/from correct Python type
  #   e.g. for Python list or NumPy or ctypes

  def scatter_atoms(self,name,type,count,data):
    if name: name = name.encode()
    self.lib.lammps_scatter_atoms(self.lmp,name,type,count,data)

  # -------------------------------------------------------------------------

  def scatter_atoms_subset(self,name,type,count,ndata,ids,data):
    if name: name = name.encode()
    self.lib.lammps_scatter_atoms_subset(self.lmp,name,type,count,ndata,ids,data)

  # return vector of atom/compute/fix properties gathered across procs
  # 3 variants to match src/library.cpp
  # name = atom property recognized by LAMMPS in atom->extract()
  # type = 0 for integer values, 1 for double values
  # count = number of per-atom valus, 1 for type or charge, 3 for x or f
  # returned data is a 1d vector - doc how it is ordered?
  # NOTE: need to insure are converting to/from correct Python type
  #   e.g. for Python list or NumPy or ctypes
  def gather(self,name,type,count):
    if name: name = name.encode()
    natoms = self.lib.lammps_get_natoms(self.lmp)
    if type == 0:
      data = ((count*natoms)*c_int)()
      self.lib.lammps_gather(self.lmp,name,type,count,data)
    elif type == 1:
      data = ((count*natoms)*c_double)()
      self.lib.lammps_gather(self.lmp,name,type,count,data)
    else: return None
    return data

  def gather_concat(self,name,type,count):
    if name: name = name.encode()
    natoms = self.lib.lammps_get_natoms(self.lmp)
    if type == 0:
      data = ((count*natoms)*c_int)()
      self.lib.lammps_gather_concat(self.lmp,name,type,count,data)
    elif type == 1:
      data = ((count*natoms)*c_double)()
      self.lib.lammps_gather_concat(self.lmp,name,type,count,data)
    else: return None
    return data

  def gather_subset(self,name,type,count,ndata,ids):
    if name: name = name.encode()
    if type == 0:
      data = ((count*ndata)*c_int)()
      self.lib.lammps_gather_subset(self.lmp,name,type,count,ndata,ids,data)
    elif type == 1:
      data = ((count*ndata)*c_double)()
      self.lib.lammps_gather_subset(self.lmp,name,type,count,ndata,ids,data)
    else: return None
    return data

  # scatter vector of atom/compute/fix properties across procs
  # 2 variants to match src/library.cpp
  # name = atom property recognized by LAMMPS in atom->extract()
  # type = 0 for integer values, 1 for double values
  # count = number of per-atom valus, 1 for type or charge, 3 for x or f
  # assume data is of correct type and length, as created by gather_atoms()
  # NOTE: need to insure are converting to/from correct Python type
  #   e.g. for Python list or NumPy or ctypes

  def scatter(self,name,type,count,data):
    if name: name = name.encode()
    self.lib.lammps_scatter(self.lmp,name,type,count,data)

  def scatter_subset(self,name,type,count,ndata,ids,data):
    if name: name = name.encode()
    self.lib.lammps_scatter_subset(self.lmp,name,type,count,ndata,ids,data)

   # -------------------------------------------------------------------------

  def encode_image_flags(self,ix,iy,iz):
    """ convert 3 integers with image flags for x-, y-, and z-direction
    into a single integer like it is used internally in LAMMPS

    This method is a wrapper around the :cpp:func:`lammps_encode_image_flags`
    function of library interface.

    :param ix: x-direction image flag
    :type  ix: int
    :param iy: y-direction image flag
    :type  iy: int
    :param iz: z-direction image flag
    :type  iz: int
    :return: encoded image flags
    :rtype: lammps.c_imageint
    """
    return self.lib.lammps_encode_image_flags(ix,iy,iz)

  # -------------------------------------------------------------------------

  def decode_image_flags(self,image):
    """ Convert encoded image flag integer into list of three regular integers.

    This method is a wrapper around the :cpp:func:`lammps_decode_image_flags`
    function of library interface.

    :param image: encoded image flags
    :type image:  lammps.c_imageint
    :return: list of three image flags in x-, y-, and z- direction
    :rtype: list of 3 int
    """

    flags = (c_int*3)()
    self.lib.lammps_decode_image_flags(image,byref(flags))

    return [int(i) for i in flags]

  # -------------------------------------------------------------------------

  # create N atoms on all procs
  # N = global number of atoms
  # id = ID of each atom (optional, can be None)
  # type = type of each atom (1 to Ntypes) (required)
  # x = coords of each atom as (N,3) array (required)
  # v = velocity of each atom as (N,3) array (optional, can be None)
  # NOTE: how could we insure are passing correct type to LAMMPS
  #   e.g. for Python list or NumPy, etc
  #   ditto for gather_atoms() above

  def create_atoms(self,n,id,type,x,v=None,image=None,shrinkexceed=False):
    """
    Create N atoms from list of coordinates and properties

    This function is a wrapper around the :cpp:func:`lammps_create_atoms`
    function of the C-library interface, and the behavior is similar except
    that the *v*, *image*, and *shrinkexceed* arguments are optional and
    default to *None*, *None*, and *False*, respectively. With none being
    equivalent to a ``NULL`` pointer in C.

    The lists of coordinates, types, atom IDs, velocities, image flags can
    be provided in any format that may be converted into the required
    internal data types.  Also the list may contain more than *N* entries,
    but not fewer.  In the latter case, the function will return without
    attempting to create atoms.  You may use the :py:func:`encode_image_flags
    <lammps.encode_image_flags>` method to properly combine three integers
    with image flags into a single integer.

    :param n: number of atoms for which data is provided
    :type n: int
    :param id: list of atom IDs with at least n elements or None
    :type id: list of lammps.tagint
    :param type: list of atom types
    :type type: list of int
    :param x: list of coordinates for x-, y-, and z (flat list of 3n entries)
    :type x: list of float
    :param v: list of velocities for x-, y-, and z (flat list of 3n entries) or None (optional)
    :type v: list of float
    :param image: list of encoded image flags (optional)
    :type image: list of lammps.imageint
    :param shrinkexceed: whether to expand shrink-wrap boundaries if atoms are outside the box (optional)
    :type shrinkexceed: bool
    :return: number of atoms created. 0 if insufficient or invalid data
    :rtype: int
    """
    if id:
      id_lmp = (self.c_tagint*n)()
      try:
        id_lmp[:] = id[0:n]
      except:
        return 0
    else:
      id_lmp = None

    type_lmp = (c_int*n)()
    try:
      type_lmp[:] = type[0:n]
    except:
      return 0

    three_n = 3*n
    x_lmp = (c_double*three_n)()
    try:
      x_lmp[:] = x[0:three_n]
    except:
      return 0

    if v:
      v_lmp = (c_double*(three_n))()
      try:
        v_lmp[:] = v[0:three_n]
      except:
        return 0
    else:
      v_lmp = None

    if image:
      img_lmp = (self.c_imageint*n)()
      try:
        img_lmp[:] = image[0:n]
      except:
        return 0
    else:
      img_lmp = None

    if shrinkexceed:
      se_lmp = 1
    else:
      se_lmp = 0

    self.lib.lammps_create_atoms.argtypes = [c_void_p, c_int, POINTER(self.c_tagint*n),
                                     POINTER(c_int*n), POINTER(c_double*three_n),
                                     POINTER(c_double*three_n),
                                     POINTER(self.c_imageint*n), c_int]
    return self.lib.lammps_create_atoms(self.lmp, n, id_lmp, type_lmp, x_lmp, v_lmp, img_lmp, se_lmp)

  # -------------------------------------------------------------------------

  @property
  def has_mpi_support(self):
    """ Report whether the LAMMPS shared library was compiled with a
    real MPI library or in serial.

    This is a wrapper around the :cpp:func:`lammps_config_has_mpi_support`
    function of the library interface.

    :return: False when compiled with MPI STUBS, otherwise True
    :rtype: bool
    """
    return self.lib.lammps_config_has_mpi_support() != 0

  # -------------------------------------------------------------------------

  @property
  def has_exceptions(self):
    """ Report whether the LAMMPS shared library was compiled with C++
    exceptions handling enabled

    This is a wrapper around the :cpp:func:`lammps_config_has_exceptions`
    function of the library interface.

    :return: state of C++ exception support
    :rtype: bool
    """
    return self.lib.lammps_config_has_exceptions() != 0

  # -------------------------------------------------------------------------

  @property
  def has_gzip_support(self):
    """ Report whether the LAMMPS shared library was compiled with support
    for reading and writing compressed files through ``gzip``.

    This is a wrapper around the :cpp:func:`lammps_config_has_gzip_support`
    function of the library interface.

    :return: state of gzip support
    :rtype: bool
    """
    return self.lib.lammps_config_has_gzip_support() != 0

  # -------------------------------------------------------------------------

  @property
  def has_png_support(self):
    """ Report whether the LAMMPS shared library was compiled with support
    for writing images in PNG format.

    This is a wrapper around the :cpp:func:`lammps_config_has_png_support`
    function of the library interface.

    :return: state of PNG support
    :rtype: bool
    """
    return self.lib.lammps_config_has_png_support() != 0

  # -------------------------------------------------------------------------

  @property
  def has_jpeg_support(self):
    """ Report whether the LAMMPS shared library was compiled with support
    for writing images in JPEG format.

    This is a wrapper around the :cpp:func:`lammps_config_has_jpeg_support`
    function of the library interface.

    :return: state of JPEG support
    :rtype: bool
    """
    return self.lib.lammps_config_has_jpeg_support() != 0

  # -------------------------------------------------------------------------

  @property
  def has_ffmpeg_support(self):
    """ State of support for writing movies with ``ffmpeg`` in the LAMMPS shared library

    This is a wrapper around the :cpp:func:`lammps_config_has_ffmpeg_support`
    function of the library interface.

    :return: state of ffmpeg support
    :rtype: bool
    """
    return self.lib.lammps_config_has_ffmpeg_support() != 0

  # -------------------------------------------------------------------------

  @property
  def installed_packages(self):
    """ List of the names of enabled packages in the LAMMPS shared library

    This is a wrapper around the functions :cpp:func:`lammps_config_package_count`
    and :cpp:func`lammps_config_package_name` of the library interface.

    :return
    """
    if self._installed_packages is None:
      self._installed_packages = []
      npackages = self.lib.lammps_config_package_count()
      sb = create_string_buffer(100)
      for idx in range(npackages):
        self.lib.lammps_config_package_name(idx, sb, 100)
        self._installed_packages.append(sb.value.decode())
    return self._installed_packages

  # -------------------------------------------------------------------------

  def has_style(self, category, name):
    """Returns whether a given style name is available in a given category

    This is a wrapper around the function :cpp:func:`lammps_has_style`
    of the library interface.

    :param category: name of category
    :type  category: string
    :param name: name of the style
    :type  name: string

    :return: true if style is available in given category
    :rtype:  bool
    """
    return self.lib.lammps_has_style(self.lmp, category.encode(), name.encode()) != 0

  # -------------------------------------------------------------------------

  def available_styles(self, category):
    """Returns a list of styles available for a given category

    This is a wrapper around the functions :cpp:func:`lammps_style_count`
    and :cpp:func`lammps_style_name` of the library interface.

    :param category: name of category
    :type  category: string

    :return: list of style names in given category
    :rtype:  list
    """
    if self._available_styles is None:
      self._available_styles = {}

    if category not in self._available_styles:
      self._available_styles[category] = []
      nstyles = self.lib.lammps_style_count(self.lmp, category.encode())
      sb = create_string_buffer(100)
      for idx in range(nstyles):
        self.lib.lammps_style_name(self.lmp, category.encode(), idx, sb, 100)
        self._available_styles[category].append(sb.value.decode())
    return self._available_styles[category]

  # -------------------------------------------------------------------------

  def set_fix_external_callback(self, fix_name, callback, caller=None):
    import numpy as np

    def _ctype_to_numpy_int(ctype_int):
          if ctype_int == c_int32:
            return np.int32
          elif ctype_int == c_int64:
            return np.int64
          return np.intc

    def callback_wrapper(caller, ntimestep, nlocal, tag_ptr, x_ptr, fext_ptr):
      tag = self.numpy.iarray(self.c_tagint, tag_ptr, nlocal, 1)
      x   = self.numpy.darray(x_ptr, nlocal, 3)
      f   = self.numpy.darray(fext_ptr, nlocal, 3)
      callback(caller, ntimestep, nlocal, tag, x, f)

    cFunc   = self.FIX_EXTERNAL_CALLBACK_FUNC(callback_wrapper)
    cCaller = caller

    self.callback[fix_name] = { 'function': cFunc, 'caller': caller }
    self.lib.lammps_set_fix_external_callback(self.lmp, fix_name.encode(), cFunc, cCaller)

  # -------------------------------------------------------------------------

  def get_neighlist(self, idx):
    """Returns an instance of :class:`NeighList` which wraps access to the neighbor list with the given index

    :param idx: index of neighbor list
    :type  idx: int
    :return: an instance of :class:`NeighList` wrapping access to neighbor list data
    :rtype:  NeighList
    """
    if idx < 0:
        return None
    return NeighList(self, idx)

  # -------------------------------------------------------------------------

  def find_pair_neighlist(self, style, exact=True, nsub=0, request=0):
    """Find neighbor list index of pair style neighbor list

    Try finding pair instance that matches style. If exact is set, the pair must
    match style exactly. If exact is 0, style must only be contained. If pair is
    of style pair/hybrid, style is instead matched the nsub-th hybrid sub-style.

    Once the pair instance has been identified, multiple neighbor list requests
    may be found. Every neighbor list is uniquely identified by its request
    index. Thus, providing this request index ensures that the correct neighbor
    list index is returned.

    :param style: name of pair style that should be searched for
    :type  style: string
    :param exact: controls whether style should match exactly or only must be contained in pair style name, defaults to True
    :type  exact: bool, optional
    :param nsub:  match nsub-th hybrid sub-style, defaults to 0
    :type  nsub:  int, optional
    :param request:   index of neighbor list request, in case there are more than one, defaults to 0
    :type  request:   int, optional
    :return: neighbor list index if found, otherwise -1
    :rtype:  int
     """
    style = style.encode()
    exact = int(exact)
    idx = self.lib.lammps_find_pair_neighlist(self.lmp, style, exact, nsub, request)
    return self.get_neighlist(idx)

  # -------------------------------------------------------------------------

  def find_fix_neighlist(self, fixid, request=0):
    """Find neighbor list index of fix neighbor list

    :param fixid: name of fix
    :type  fixid: string
    :param request:   index of neighbor list request, in case there are more than one, defaults to 0
    :type  request:   int, optional
    :return: neighbor list index if found, otherwise -1
    :rtype:  int
     """
    fixid = fixid.encode()
    idx = self.lib.lammps_find_fix_neighlist(self.lmp, fixid, request)
    return self.get_neighlist(idx)

  # -------------------------------------------------------------------------

  def find_compute_neighlist(self, computeid, request=0):
    """Find neighbor list index of compute neighbor list

    :param computeid: name of compute
    :type  computeid: string
    :param request:   index of neighbor list request, in case there are more than one, defaults to 0
    :type  request:   int, optional
    :return: neighbor list index if found, otherwise -1
    :rtype:  int
     """
    computeid = computeid.encode()
    idx = self.lib.lammps_find_compute_neighlist(self.lmp, computeid, request)
    return self.get_neighlist(idx)

  # -------------------------------------------------------------------------

  def get_neighlist_size(self, idx):
    """Return the number of elements in neighbor list with the given index

    :param idx: neighbor list index
    :type  idx: int
    :return: number of elements in neighbor list with index idx
    :rtype:  int
     """
    return self.lib.lammps_neighlist_num_elements(self.lmp, idx)

  # -------------------------------------------------------------------------

  def get_neighlist_element_neighbors(self, idx, element):
    """Return data of neighbor list entry

    :param element: neighbor list index
    :type  element: int
    :param element: neighbor list element index
    :type  element: int
    :return: tuple with atom local index, number of neighbors and array of neighbor local atom indices
    :rtype:  (int, int, numpy.array)
    """
    c_iatom = c_int()
    c_numneigh = c_int()
    c_neighbors = POINTER(c_int)()
    self.lib.lammps_neighlist_element_neighbors(self.lmp, idx, element, byref(c_iatom), byref(c_numneigh), byref(c_neighbors))
    neighbors = self.numpy.iarray(c_int, c_neighbors, c_numneigh.value, 1)
    return c_iatom.value, c_numneigh.value, neighbors

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

################################################################################
# Alternative Python Wrapper
# Written by Richard Berger <richard.berger@temple.edu>
################################################################################

class OutputCapture(object):
  """ Utility class to capture LAMMPS library output """

  def __init__(self):
    self.stdout_pipe_read, self.stdout_pipe_write = os.pipe()
    self.stdout_fd = 1

  def __enter__(self):
    self.stdout = os.dup(self.stdout_fd)
    os.dup2(self.stdout_pipe_write, self.stdout_fd)
    return self

  def __exit__(self, type, value, tracebac):
    os.dup2(self.stdout, self.stdout_fd)
    os.close(self.stdout)
    os.close(self.stdout_pipe_read)
    os.close(self.stdout_pipe_write)

  # check if we have more to read from the pipe
  def more_data(self, pipe):
    r, _, _ = select.select([pipe], [], [], 0)
    return bool(r)

  # read the whole pipe
  def read_pipe(self, pipe):
    out = ""
    while self.more_data(pipe):
      out += os.read(pipe, 1024).decode()
    return out

  @property
  def output(self):
    return self.read_pipe(self.stdout_pipe_read)

# -------------------------------------------------------------------------

class Variable(object):
  def __init__(self, lammps_wrapper_instance, name, style, definition):
    self.wrapper = lammps_wrapper_instance
    self.name = name
    self.style = style
    self.definition = definition.split()

  @property
  def value(self):
    if self.style == 'atom':
      return list(self.wrapper.lmp.extract_variable(self.name, "all", 1))
    else:
      value = self.wrapper.lmp_print('"${%s}"' % self.name).strip()
      try:
        return float(value)
      except ValueError:
        return value

# -------------------------------------------------------------------------

class AtomList(object):
  def __init__(self, lammps_wrapper_instance):
    self.lmp = lammps_wrapper_instance
    self.natoms = self.lmp.system.natoms
    self.dimensions = self.lmp.system.dimensions

  def __getitem__(self, index):
    if self.dimensions == 2:
        return Atom2D(self.lmp, index + 1)
    return Atom(self.lmp, index + 1)

# -------------------------------------------------------------------------

class Atom(object):
  def __init__(self, lammps_wrapper_instance, index):
    self.lmp = lammps_wrapper_instance
    self.index = index

  @property
  def id(self):
    return int(self.lmp.eval("id[%d]" % self.index))

  @property
  def type(self):
    return int(self.lmp.eval("type[%d]" % self.index))

  @property
  def mol(self):
    return self.lmp.eval("mol[%d]" % self.index)

  @property
  def mass(self):
    return self.lmp.eval("mass[%d]" % self.index)

  @property
  def position(self):
    return (self.lmp.eval("x[%d]" % self.index),
            self.lmp.eval("y[%d]" % self.index),
            self.lmp.eval("z[%d]" % self.index))

  @position.setter
  def position(self, value):
     self.lmp.set("atom", self.index, "x", value[0])
     self.lmp.set("atom", self.index, "y", value[1])
     self.lmp.set("atom", self.index, "z", value[2])

  @property
  def velocity(self):
    return (self.lmp.eval("vx[%d]" % self.index),
            self.lmp.eval("vy[%d]" % self.index),
            self.lmp.eval("vz[%d]" % self.index))

  @velocity.setter
  def velocity(self, value):
     self.lmp.set("atom", self.index, "vx", value[0])
     self.lmp.set("atom", self.index, "vy", value[1])
     self.lmp.set("atom", self.index, "vz", value[2])

  @property
  def force(self):
    return (self.lmp.eval("fx[%d]" % self.index),
            self.lmp.eval("fy[%d]" % self.index),
            self.lmp.eval("fz[%d]" % self.index))

  @property
  def charge(self):
    return self.lmp.eval("q[%d]" % self.index)

# -------------------------------------------------------------------------

class Atom2D(Atom):
  def __init__(self, lammps_wrapper_instance, index):
    super(Atom2D, self).__init__(lammps_wrapper_instance, index)

  @property
  def position(self):
    return (self.lmp.eval("x[%d]" % self.index),
            self.lmp.eval("y[%d]" % self.index))

  @position.setter
  def position(self, value):
     self.lmp.set("atom", self.index, "x", value[0])
     self.lmp.set("atom", self.index, "y", value[1])

  @property
  def velocity(self):
    return (self.lmp.eval("vx[%d]" % self.index),
            self.lmp.eval("vy[%d]" % self.index))

  @velocity.setter
  def velocity(self, value):
     self.lmp.set("atom", self.index, "vx", value[0])
     self.lmp.set("atom", self.index, "vy", value[1])

  @property
  def force(self):
    return (self.lmp.eval("fx[%d]" % self.index),
            self.lmp.eval("fy[%d]" % self.index))

# -------------------------------------------------------------------------

class variable_set:
    def __init__(self, name, variable_dict):
        self._name = name
        array_pattern = re.compile(r"(?P<arr>.+)\[(?P<index>[0-9]+)\]")

        for key, value in variable_dict.items():
            m = array_pattern.match(key)
            if m:
                g = m.groupdict()
                varname = g['arr']
                idx = int(g['index'])
                if varname not in self.__dict__:
                    self.__dict__[varname] = {}
                self.__dict__[varname][idx] = value
            else:
                self.__dict__[key] = value

    def __str__(self):
        return "{}({})".format(self._name, ','.join(["{}={}".format(k, self.__dict__[k]) for k in self.__dict__.keys() if not k.startswith('_')]))

    def __repr__(self):
        return self.__str__()

# -------------------------------------------------------------------------

def get_thermo_data(output):
    """ traverse output of runs and extract thermo data columns """
    if isinstance(output, str):
        lines = output.splitlines()
    else:
        lines = output

    runs = []
    columns = []
    in_run = False
    current_run = {}

    for line in lines:
        if line.startswith("Per MPI rank memory allocation"):
            in_run = True
        elif in_run and len(columns) == 0:
            # first line after memory usage are column names
            columns = line.split()

            current_run = {}

            for col in columns:
                current_run[col] = []

        elif line.startswith("Loop time of "):
            in_run = False
            columns = None
            thermo_data = variable_set('ThermoData', current_run)
            r = {'thermo' : thermo_data }
            runs.append(namedtuple('Run', list(r.keys()))(*list(r.values())))
        elif in_run and len(columns) > 0:
            items = line.split()
            # Convert thermo output and store it.
            # It must have the same number of columns and
            # all of them must be convertible to floats.
            # Otherwise we ignore the line
            if len(items) == len(columns):
                try:
                    values = [float(x) for x in items]
                    for i, col in enumerate(columns):
                        current_run[col].append(values[i])
                except ValueError:
                  pass

    return runs

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

class PyLammps(object):
  """
  More Python-like wrapper for LAMMPS (e.g., for IPython)
  See examples/ipython for usage
  """

  def __init__(self,name="",cmdargs=None,ptr=None,comm=None):
    if ptr:
      if isinstance(ptr,PyLammps):
        self.lmp = ptr.lmp
      elif isinstance(ptr,lammps):
        self.lmp = ptr
      else:
        self.lmp = lammps(name=name,cmdargs=cmdargs,ptr=ptr,comm=comm)
    else:
      self.lmp = lammps(name=name,cmdargs=cmdargs,ptr=None,comm=comm)
    print("LAMMPS output is captured by PyLammps wrapper")
    self._cmd_history = []
    self.runs = []

  def __del__(self):
    if self.lmp: self.lmp.close()
    self.lmp = None

  def close(self):
    if self.lmp: self.lmp.close()
    self.lmp = None

  def version(self):
    return self.lmp.version()

  def file(self,file):
    self.lmp.file(file)

  def write_script(self,filename):
    """ Write LAMMPS script file containing all commands executed up until now """
    with open(filename, "w") as f:
      for cmd in self._cmd_history:
        f.write("%s\n" % cmd)

  def command(self,cmd):
    self.lmp.command(cmd)
    self._cmd_history.append(cmd)

  def run(self, *args, **kwargs):
    output = self.__getattr__('run')(*args, **kwargs)

    if(self.has_mpi4py):
      output = self.lmp.comm.bcast(output, root=0)

    self.runs += get_thermo_data(output)
    return output

  @property
  def last_run(self):
    if len(self.runs) > 0:
        return self.runs[-1]
    return None

  @property
  def atoms(self):
    return AtomList(self)

  @property
  def system(self):
    output = self.info("system")
    d = self._parse_info_system(output)
    return namedtuple('System', d.keys())(*d.values())

  @property
  def communication(self):
    output = self.info("communication")
    d = self._parse_info_communication(output)
    return namedtuple('Communication', d.keys())(*d.values())

  @property
  def computes(self):
    output = self.info("computes")
    return self._parse_element_list(output)

  @property
  def dumps(self):
    output = self.info("dumps")
    return self._parse_element_list(output)

  @property
  def fixes(self):
    output = self.info("fixes")
    return self._parse_element_list(output)

  @property
  def groups(self):
    output = self.info("groups")
    return self._parse_groups(output)

  @property
  def variables(self):
    output = self.info("variables")
    vars = {}
    for v in self._parse_element_list(output):
      vars[v['name']] = Variable(self, v['name'], v['style'], v['def'])
    return vars

  def eval(self, expr):
    value = self.lmp_print('"$(%s)"' % expr).strip()
    try:
      return float(value)
    except ValueError:
      return value

  def _split_values(self, line):
    return [x.strip() for x in line.split(',')]

  def _get_pair(self, value):
    return [x.strip() for x in value.split('=')]

  def _parse_info_system(self, output):
    lines = output[6:-2]
    system = {}

    for line in lines:
      if line.startswith("Units"):
        system['units'] = self._get_pair(line)[1]
      elif line.startswith("Atom style"):
        system['atom_style'] = self._get_pair(line)[1]
      elif line.startswith("Atom map"):
        system['atom_map'] = self._get_pair(line)[1]
      elif line.startswith("Atoms"):
        parts = self._split_values(line)
        system['natoms'] = int(self._get_pair(parts[0])[1])
        system['ntypes'] = int(self._get_pair(parts[1])[1])
        system['style'] = self._get_pair(parts[2])[1]
      elif line.startswith("Kspace style"):
        system['kspace_style'] = self._get_pair(line)[1]
      elif line.startswith("Dimensions"):
        system['dimensions'] = int(self._get_pair(line)[1])
      elif line.startswith("Orthogonal box"):
        system['orthogonal_box'] = [float(x) for x in self._get_pair(line)[1].split('x')]
      elif line.startswith("Boundaries"):
        system['boundaries'] = self._get_pair(line)[1]
      elif line.startswith("xlo"):
        keys, values = [self._split_values(x) for x in self._get_pair(line)]
        for key, value in zip(keys, values):
          system[key] = float(value)
      elif line.startswith("ylo"):
        keys, values = [self._split_values(x) for x in self._get_pair(line)]
        for key, value in zip(keys, values):
          system[key] = float(value)
      elif line.startswith("zlo"):
        keys, values = [self._split_values(x) for x in self._get_pair(line)]
        for key, value in zip(keys, values):
          system[key] = float(value)
      elif line.startswith("Molecule type"):
        system['molecule_type'] = self._get_pair(line)[1]
      elif line.startswith("Bonds"):
        parts = self._split_values(line)
        system['nbonds'] = int(self._get_pair(parts[0])[1])
        system['nbondtypes'] = int(self._get_pair(parts[1])[1])
        system['bond_style'] = self._get_pair(parts[2])[1]
      elif line.startswith("Angles"):
        parts = self._split_values(line)
        system['nangles'] = int(self._get_pair(parts[0])[1])
        system['nangletypes'] = int(self._get_pair(parts[1])[1])
        system['angle_style'] = self._get_pair(parts[2])[1]
      elif line.startswith("Dihedrals"):
        parts = self._split_values(line)
        system['ndihedrals'] = int(self._get_pair(parts[0])[1])
        system['ndihedraltypes'] = int(self._get_pair(parts[1])[1])
        system['dihedral_style'] = self._get_pair(parts[2])[1]
      elif line.startswith("Impropers"):
        parts = self._split_values(line)
        system['nimpropers'] = int(self._get_pair(parts[0])[1])
        system['nimpropertypes'] = int(self._get_pair(parts[1])[1])
        system['improper_style'] = self._get_pair(parts[2])[1]

    return system

  def _parse_info_communication(self, output):
    lines = output[6:-3]
    comm = {}

    for line in lines:
      if line.startswith("MPI library"):
        comm['mpi_version'] = line.split(':')[1].strip()
      elif line.startswith("Comm style"):
        parts = self._split_values(line)
        comm['comm_style'] = self._get_pair(parts[0])[1]
        comm['comm_layout'] = self._get_pair(parts[1])[1]
      elif line.startswith("Processor grid"):
        comm['proc_grid'] = [int(x) for x in self._get_pair(line)[1].split('x')]
      elif line.startswith("Communicate velocities for ghost atoms"):
        comm['ghost_velocity'] = (self._get_pair(line)[1] == "yes")
      elif line.startswith("Nprocs"):
        parts = self._split_values(line)
        comm['nprocs'] = int(self._get_pair(parts[0])[1])
        comm['nthreads'] = int(self._get_pair(parts[1])[1])
    return comm

  def _parse_element_list(self, output):
    lines = output[6:-3]
    elements = []

    for line in lines:
      element_info = self._split_values(line.split(':')[1].strip())
      element = {'name': element_info[0]}
      for key, value in [self._get_pair(x) for x in element_info[1:]]:
        element[key] = value
      elements.append(element)
    return elements

  def _parse_groups(self, output):
    lines = output[6:-3]
    groups = []
    group_pattern = re.compile(r"(?P<name>.+) \((?P<type>.+)\)")

    for line in lines:
      m = group_pattern.match(line.split(':')[1].strip())
      group = {'name': m.group('name'), 'type': m.group('type')}
      groups.append(group)
    return groups

  def lmp_print(self, s):
    """ needed for Python2 compatibility, since print is a reserved keyword """
    return self.__getattr__("print")(s)

  def __dir__(self):
    return ['angle_coeff', 'angle_style', 'atom_modify', 'atom_style', 'atom_style',
    'bond_coeff', 'bond_style', 'boundary', 'change_box', 'communicate', 'compute',
    'create_atoms', 'create_box', 'delete_atoms', 'delete_bonds', 'dielectric',
    'dihedral_coeff', 'dihedral_style', 'dimension', 'dump', 'fix', 'fix_modify',
    'group', 'improper_coeff', 'improper_style', 'include', 'kspace_modify',
    'kspace_style', 'lattice', 'mass', 'minimize', 'min_style', 'neighbor',
    'neigh_modify', 'newton', 'nthreads', 'pair_coeff', 'pair_modify',
    'pair_style', 'processors', 'read', 'read_data', 'read_restart', 'region',
    'replicate', 'reset_timestep', 'restart', 'run', 'run_style', 'thermo',
    'thermo_modify', 'thermo_style', 'timestep', 'undump', 'unfix', 'units',
    'variable', 'velocity', 'write_restart']

  def __getattr__(self, name):
    def handler(*args, **kwargs):
      cmd_args = [name] + [str(x) for x in args]

      with OutputCapture() as capture:
        self.command(' '.join(cmd_args))
        output = capture.output

      if 'verbose' in kwargs and kwargs['verbose']:
        print(output)

      lines = output.splitlines()

      if len(lines) > 1:
        return lines
      elif len(lines) == 1:
        return lines[0]
      return None

    return handler


class IPyLammps(PyLammps):
  """
  IPython wrapper for LAMMPS which adds embedded graphics capabilities
  """

  def __init__(self,name="",cmdargs=None,ptr=None,comm=None):
    super(IPyLammps, self).__init__(name=name,cmdargs=cmdargs,ptr=ptr,comm=comm)

  def image(self, filename="snapshot.png", group="all", color="type", diameter="type",
            size=None, view=None, center=None, up=None, zoom=1.0):
    cmd_args = [group, "image", filename, color, diameter]

    if size:
      width = size[0]
      height = size[1]
      cmd_args += ["size", width, height]

    if view:
      theta = view[0]
      phi = view[1]
      cmd_args += ["view", theta, phi]

    if center:
      flag = center[0]
      Cx = center[1]
      Cy = center[2]
      Cz = center[3]
      cmd_args += ["center", flag, Cx, Cy, Cz]

    if up:
      Ux = up[0]
      Uy = up[1]
      Uz = up[2]
      cmd_args += ["up", Ux, Uy, Uz]

    if zoom:
      cmd_args += ["zoom", zoom]

    cmd_args.append("modify backcolor white")

    self.write_dump(*cmd_args)
    from IPython.core.display import Image
    return Image('snapshot.png')

  def video(self, filename):
    from IPython.display import HTML
    return HTML("<video controls><source src=\"" + filename + "\"></video>")
