.. index:: angle\_style cosine/delta

angle\_style cosine/delta command
=================================

angle\_style cosine/delta/omp command
=====================================

Syntax
""""""


.. parsed-literal::

   angle_style cosine/delta

Examples
""""""""


.. parsed-literal::

   angle_style cosine/delta
   angle_coeff 2\*4 75.0 100.0

Description
"""""""""""

The *cosine/delta* angle style uses the potential

.. image:: Eqs/angle_cosine_delta.jpg
   :align: center

where theta0 is the equilibrium value of the angle, and K is a
prefactor.  Note that the usual 1/2 factor is included in K.

The following coefficients must be defined for each angle type via the
:doc:`angle\_coeff <angle_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read\_data <read_data>`
or :doc:`read\_restart <read_restart>` commands:

* K (energy)
* theta0 (degrees)

Theta0 is specified in degrees, but LAMMPS converts it to radians
internally.


----------


Styles with a *gpu*\ , *intel*\ , *kk*\ , *omp*\ , or *opt* suffix are
functionally the same as the corresponding style without the suffix.
They have been optimized to run faster, depending on your available
hardware, as discussed on the :doc:`Speed packages <Speed_packages>` doc
page.  The accelerated styles take the same arguments and should
produce the same results, except for round-off and precision issues.

These accelerated styles are part of the GPU, USER-INTEL, KOKKOS,
USER-OMP and OPT packages, respectively.  They are only enabled if
LAMMPS was built with those packages.  See the :doc:`Build package <Build_package>` doc page for more info.

You can specify the accelerated styles explicitly in your input script
by including their suffix, or you can use the :doc:`-suffix command-line switch <Run_options>` when you invoke LAMMPS, or you can use the
:doc:`suffix <suffix>` command in your input script.

See the :doc:`Speed packages <Speed_packages>` doc page for more
instructions on how to use the accelerated styles effectively.


----------


Restrictions
""""""""""""


This angle style can only be used if LAMMPS was built with the
MOLECULE package.  See the :doc:`Build package <Build_package>` doc page
for more info.

Related commands
""""""""""""""""

:doc:`angle\_coeff <angle_coeff>`, :doc:`angle\_style cosine/squared <angle_cosine_squared>`

**Default:** none


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
