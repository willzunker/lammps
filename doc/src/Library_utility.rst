Library interface utility functions
===================================

To simplify some of the tasks, the library interface contains
some utility functions that are not directly calling LAMMPS.

-----------------------

.. doxygenfunction:: lammps_encode_image_flags
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_decode_image_flags(int image, int *flags)
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_set_fix_external_callback(void *, char *, FixExternalFnPtr, void*)
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_has_error
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_get_last_error_message
   :project: progguide
