.. index:: fix numdiff

fix numdiff command
====================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID numdiff Nevery Delta

* ID, group-ID are documented in :doc:`fix <fix>` command
* numdiff = style name of this fix command
* Nevery = calculate force by finite difference every this many timesteps
* Delta = finite difference displacement length (distance units)
  

Examples
""""""""


.. parsed-literal::

   fix 1 all numdiff 1 0.0001
   fix 1 all numdiff 10 1e-6
   fix 1 all numdiff 100 0.01

Description
"""""""""""

Calculate forces through finite difference calculations of energy
versus position.  These forces can be compared to analytic forces
computed by pair styles, bond styles, etc.  E.g. for debugging
purposes.

The group specified with the command means only atoms within the group
have their averages computed.  Results are set to 0.0 for atoms not in
the group.

This fix performs a loop over all atoms (in the group).  For each atom
and each component of force it adds *Delta* to the position, and
computes the new energy of the entire system.  It then subtracts
*Delta* from the original position and again computes the new energy
of the system.  It then restores the original position.  That
component of force is calculated as the difference in energy divided
by two times *Delta*.

.. note::

   It is important to choose a suitable value for delta, the magnitude of
   atom displacements that are used to generate finite difference
   approximations to the exact forces.  For typical systems, a value in
   the range of 1 part in 1e4 to 1e5 of the typical separation distance
   between atoms in the liquid or solid state will be sufficient.
   However, the best value will depend on a multitude of factors
   including the stiffness of the interatomic potential, the thermodynamic
   state of the material being probed, and so on. The only way to be sure
   that you have made a good choice is to do a sensitivity study on a
   representative atomic configuration, sweeping over a wide range of
   values of delta.  If delta is too small, the output forces will vary
   erratically due to truncation effects. If delta is increased beyond a
   certain point, the output forces will start to vary smoothly with
   delta, due to growing contributions from higher order derivatives. In
   between these two limits, the numerical force values should be largely
   independent of delta.

.. note::

   The cost of each energy evaluation is essentially the cost of an MD
   timestep.  This invoking this fix once has a cost of 2N timesteps,
   where N is the total number of atoms in the system (assuming all atoms
   are included in the group).  So this fix can be very expensive to use
   for large systems.

----------


The *Nevery*\ argument specifies on what timesteps the force will 
be used calculated by finite difference.

The *Delta*\ argument specifies the positional displacement each
atom will undergo.


----------


**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files
<restart>`.  None of the :doc:`fix_modify <fix_modify>` options are
relevant to this fix.

This fix produces a per-atom array which can be accessed by various
:doc:`output commands <Howto_output>`, which stores the components of
the force on each atom as calculated by finite difference.  The
per-atom values can only be accessed on timesteps that are multiples
of *Nfreq* since that is when the finite difference forces are
calculated.

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is invoked during :doc:`energy
minimization <minimize>`.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`dynamical_matrix <dynamical_matrix>`,

**Default:** none
