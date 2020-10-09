Utility functions
=================

To simplify some tasks, the library interface contains these utility
functions.  They do not directly call the LAMMPS library.

- :cpp:func:`lammps_encode_image_flags`
- :cpp:func:`lammps_decode_image_flags`
- :cpp:func:`lammps_set_fix_external_callback`
- :cpp:func:`lammps_fix_external_set_energy_global`
- :cpp:func:`lammps_fix_external_set_virial_global`
- :cpp:func:`lammps_is_running`
- :cpp:func:`lammps_force_timeout`
- :cpp:func:`lammps_has_error`
- :cpp:func:`lammps_get_last_error_message`

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

.. doxygenfunction:: lammps_fix_external_set_energy_global
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_fix_external_set_virial_global
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_is_running
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_force_timeout
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_has_error
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_get_last_error_message
   :project: progguide
