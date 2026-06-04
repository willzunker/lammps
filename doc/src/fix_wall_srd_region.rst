.. index:: fix wall/srd/region

fix wall/srd/region command
===========================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID wall/srd/region region-ID keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* wall/srd/region = style name of this fix command
* region-ID = ID of a region the SRD particles should be confined to
* zero or more keyword/value pairs may be appended
* keyword = *phantom*

  .. parsed-literal::

       *phantom* value = *yes* or *no*
         *yes* = inject Lamura/Gompper virtual ("phantom") particles in
                 boundary cells (default; required for correct no-slip
                 at curved walls in the dense MPCD regime)
         *no*  = bare bounce-back only (for diagnostic comparison)

Examples
""""""""

.. code-block:: LAMMPS

   # Sphere confining the gas inside it
   region cavity sphere 0 0 0 5 side in
   fix walls all wall/srd/region cavity

   # Cylinder with shrinking radius (piston-style compression)
   variable R equal vdisplace(5.0,-0.1)
   region pipe cylinder z 0 0 v_R -8 8 side in
   fix walls all wall/srd/region pipe phantom yes

   # Bi-convex tablet die (cylinder intersected with two cup spheres)
   region die         cylinder z 0 0 4e-3 -5e-3 2e-2 side in units box
   region lowerSphere sphere 0 0 v_lower_z 6.67e-3 side in units box
   region upperSphere sphere 0 0 v_upper_z 6.67e-3 side in units box
   region cavity      intersect 3 die lowerSphere upperSphere
   fix walls all wall/srd/region cavity phantom yes

Description
"""""""""""

Confine SRD (stochastic rotation dynamics) particles inside (or outside)
an arbitrary LAMMPS :doc:`region <region>`.  This is the analogue of
:doc:`fix wall/srd <fix_wall_srd>` for non-planar geometries: instead of
specifying up to six axis-aligned faces, the wall geometry is delegated
to the Region object.  Any region style works -- ``sphere``,
``cylinder``, ``cone``, ``ellipsoid``, and compositional ``intersect``
and ``union`` regions -- as does any motion (``move``, ``rotate``) or
time-varying shape (variable radius) the region supports.  This mirrors
the pattern of :doc:`fix wall/gran/region <fix_wall_gran_region>` for
granular contact.

The wall interaction is invoked by the :doc:`fix srd <fix_srd>` command,
not by this fix directly: only the group of SRD particles tracked by fix
srd is affected, and the group-ID argument here is ignored.

A particle/wall collision occurs whenever an SRD particle would move
across the region surface in a time step.  The fix detects this by
testing whether the SRD has crossed (i.e. is no longer in the region's
``match`` half-space) and, if so, pushes the SRD back to the surface at
the end of the step, then reflects its velocity according to the
*collision* style chosen on fix srd (``slip`` or ``noslip``).

The reflection velocity calculation accounts for:

* The wall's translational velocity (region ``move`` keyword)
* The wall's rotational velocity (region ``rotate`` keyword)
* The radial expansion/contraction rate for variable-shape regions (e.g.
  a sphere or cylinder with a variable radius)

This is done by querying ``Region::velocity_contact()`` at the contact
point each step.

The region can be either *interior* (``side in``, the default, gas
inside the region) or *exterior* (``side out``, gas outside the region).
For a bi-convex tablet die the cup spheres use ``side in`` (gas inside
the sphere, sphere center in the gas cavity).  For a bi-concave shape
(convex bump punches), the bump spheres use ``side out`` (gas outside
the sphere, sphere center in the punch material).

Multiple instances of this fix can be used to confine the SRD gas with
multiple primitive regions simultaneously; alternatively a single
compound (``intersect`` / ``union``) region can be used with one fix.

Phantom particles
^^^^^^^^^^^^^^^^^

The *phantom* keyword controls whether the Lamura/Gompper virtual
("phantom") particle fill is performed in cells that the wall slices.
Phantom particles are needed at **curved** no-slip walls to avoid the
well-known density and temperature artifacts of naive bounce-back in the
dense MPCD regime (Lamura, Gompper, Ihle, Kroll, EPL 56:319, 2001;
Bolintineanu, Lechman, Plimpton, Grest, PRE 86:066703, 2012 -- the
latter cited in :doc:`fix srd <fix_srd>`).

For each SRD collision cell that the wall slices, the fix injects enough
virtual particles into the "outside-the-fluid" volume of the cell to
restore the total cell occupancy to the bulk value.  Each virtual
particle's velocity is sampled from a Maxwell-Boltzmann distribution at
the wall's Tsrd, with a mean equal to the wall's velocity at the
virtual's position.  The virtuals participate in that cell's velocity
rotation, then are discarded -- they are never streamed.

The bulk target density is computed automatically at setup time as
``NSRD * V_cell / V_fluid`` where ``V_fluid`` is the volume of the
confining region (estimated via Monte Carlo sampling).  This differs
from fix srd's ``srd_per_cell`` (which divides NSRD by all cells in the
simulation box, including wall-material cells) and gives the correct
bulk density to enforce.

For planar walls aligned with the simulation box axes the Ihle-Kroll
random grid shift (``shift yes`` on fix srd) already eliminates the
density artifact; phantom particles are not needed there and so
:doc:`fix wall/srd <fix_wall_srd>` does not implement them.  For curved
walls the shift averages the artifact but does not eliminate it, so
phantom particles are required for correct no-slip.  The default is
*yes*; set to *no* only for diagnostic comparison.

Output info
"""""""""""

This fix computes a global vector of length 3 containing the x, y, z
components of the net force on the wall from the SRDs over the most
recent collision step.  It can be accessed by various :doc:`output
commands <Howto_output>` as ``f_id[1]``, ``f_id[2]``, ``f_id[3]``.  The
vector values are "extensive".

The net force on a static, symmetric wall (e.g. a centered cylinder or
sphere with isotropic gas) is zero by symmetry; the values reported here
are then the per-step stochastic momentum-transfer vector.  To recover
the magnitude of the force (and hence the gas pressure on the wall),
time-average with :doc:`fix ave/time <fix_ave_time>`.

For a moving wall, the time-averaged force is the gas-on-wall drag in
the direction opposing motion.  For a compressing wall (varshape region
with shrinking dimension), the time-averaged force grows with the gas
pressure as the cavity volume shrinks.

No information about this fix is written to :doc:`binary restart files
<restart>`.  None of the :doc:`fix_modify <fix_modify>` options are
relevant to this fix.  No parameter of this fix can be used with the
*start/stop* keywords of the :doc:`run <run>` command.  This fix is not
invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix is part of the SRD package; see the :doc:`Build package
<Build_package>` doc page.

This fix can only be used with fix srd.  It must be defined after the
fix srd command -- LAMMPS will issue an error otherwise.

Triclinic simulation boxes are not yet supported.

The compound regions (``intersect`` / ``union``) work; the underlying
fix uses the region's surface query to find a push-back contact when the
SRD has escaped.  At the rim where two sub-region boundaries meet, the
compound's ``surface_exterior`` may return no contact (it filters by the
other sub-regions' match test); when that happens the fix falls back to
walking the sub-regions manually and picks the closest push-back from
any sub-region the SRD has actually escaped.  A defensive cap of 50
bounces per SRD per step prevents pathological geometries from
infinite-looping.

There is currently no exact-time collision solver for region walls; the
inexact "push-to-surface at end of step" reflection is always used.
When combined with ``fix srd overlap yes``, the relative ordering of
region-wall collisions vs. sphere/wall exact-time collisions is
therefore approximate, although each wall reflection itself is correct.
For slow wall motion (typical of quasi-static compaction), this is
inconsequential.

Related commands
""""""""""""""""

:doc:`fix srd <fix_srd>`, :doc:`fix wall/srd <fix_wall_srd>`,
:doc:`fix wall/gran/region <fix_wall_gran_region>`,
:doc:`region <region>`

Default
"""""""

The default for the *phantom* keyword is *yes*.
