LAMMPS Library Interfaces
*************************

As described on the :doc:`library interface to LAMMPS <Howto_library>`
doc page, LAMMPS can be built as a library (static or shared), so that
it can be called by another code, used in a :doc:`coupled manner
<Howto_couple>` with other codes, or driven through a :doc:`Python
script <Python_head>`.  Even the LAMMPS standalone executable is
essentially a thin wrapper on top of the LAMMPS library, creating a
LAMMPS instance, processing input and then existing.

Several of these approaches are based on C language wrapper functions
in the files ``src/library.h`` and ``src/library.cpp``, but it is also
possible to use C++ directly.  The basic procedure is always the same:
you create one or more instances of the
:cpp:class:`LAMMPS <LAMMPS_NS::LAMMPS>` and then pass commands as
strings or from files to that LAMMPS instance to execute calculations,
or read, manipulate, and update data from the active class instances
inside the LAMMPS to do analysis or perform operations that are not
possible with existing commands.

.. _thread-safety:

.. admonition:: Thread-safety
   :class: note

   LAMMPS was initially not conceived as a thread-safe program, but over
   the years changes have been applied to replace operations that
   collide with creating multiple LAMMPS instances from multiple-threads
   of the same process with thread-safe alternatives.  This primarily
   applies to the core LAMMPS code and less so on add-on packages,
   especially when those packages require additional code in the *lib*
   folder, interface LAMMPS to Fortran libraries, or the code uses
   static variables (like the USER-COLVARS package).

   Another major issue to deal with is to correctly handle MPI.
   Creating a LAMMPS instance requires passing an MPI communicator, or
   it assumes the ``MPI_COMM_WORLD`` communicator, which spans all MPI
   processor ranks.  When creating multiple LAMMPS object instances from
   different threads, this communicator has to be different for each
   thread or else collisions can happen.  or it has to be guaranteed,
   that only one thread at a time is active.  MPI communicators,
   however, are not a problem, if LAMMPS is compiled with the MPI STUBS
   library, which implies that there is no MPI communication and only 1
   MPI rank.

----------

.. _lammps_c_api:

LAMMPS C Library API
====================

The C library interface is most commonly used path to manage LAMMPS
instances from a compiled code and it is the basis for the :doc:`Python
<pg_python>` and :doc:`Fortran <pg_fortran>` modules.  Almost all
functions of the C language API require an argument containing a
"handle" in the form of a ``void *`` type variable, which points to the
location of a LAMMPS class instance.

The ``library.h`` header file by default includes the ``mpi.h`` header
for an MPI library, so it must be present when compiling code using the
library interface.  This usually must be the header from the same MPI
library as the LAMMPS library was compiled with.  The exception is when
LAMMPS was compiled in serial mode using the ``STUBS`` MPI library.  In
that case the calling code may be compiled with a different MPI library
for as long as :cpp:func:`lammps_open_no_mpi` is called to create a
LAMMPS instance. Then you may set the define ``-DLAMMPS_LIB_NO_MPI``
when compiling your code and the inclusion of ``mpi.h`` will be skipped
and consequently the function :cpp:func:`lammps_open` may not be used.

.. admonition:: Errors versus exceptions
   :class: note

   If any of the function calls in the LAMMPS library API will trigger
   an error inside LAMMPS, this will result in an abort of the entire
   program.  This is not always desirable.  Instead, LAMMPS can be
   compiled to instead :ref:`throw a C++ exception <exceptions>`.

.. warning::

   No checks are made on the arguments of the function calls of the C
   library interface.  *All* function arguments must be non-NULL unless
   *explicitly* allowed and point to consistent and valid data.  Buffers
   for storing returned data must be allocated to a suitable size.
   Passing invalid or unsuitable information will likely cause crashes
   or corrupt data.

------------------------------

.. toctree::
   :maxdepth: 1

   pg_lib_create
   pg_lib_execute
   pg_lib_properties
   pg_lib_objects
   pg_lib_scatter
   pg_lib_neighbor
   pg_lib_config
   pg_lib_utility
   pg_lib_add

--------------------

.. _lammps_python_api:

LAMMPS Python APIs
==================

The LAMMPS Python module enables calling the LAMMPS C library API from
Python by dynamically loading functions in the LAMMPS shared library through
the `Python ctypes module <https://docs.python.org/3/library/ctypes.html>`_.
Because of the dynamic loading, it is **required** that LAMMPS is compiled
in :ref:`"shared" mode <exe>`.  The Python interface is object oriented, but
otherwise trying to be very similar to the C library API.  Three different
Python classes to run LAMMPS are available and they build on each other.

.. toctree::
   :maxdepth: 1

   pg_python

-------------------

.. _lammps_fortran_api:

LAMMPS Fortran API
==================

The LAMMPS Fortran module is a wrapper around calling functions from the
LAMMPS C library API from Fortran through the ISO_C_BINDING feature in
Fortran 2003.  The interface is object oriented but otherwise trying to
be very similar to the C library API and the basic Python module.

.. toctree::
   :maxdepth: 1

   pg_fortran

-------------------

.. _lammps_cplusplus_api:

LAMMPS C++ API
==============

It is also possible to invoke the LAMMPS C++ API directly in your code.
It is lacking some of the convenience of the C library API, but it allows
a more direct access to simulation data and thus more low-level manipulations.
The following links provide some examples and references to the C++ API.

.. toctree::
   :maxdepth: 1

   pg_cplusplus


