.. index:: pair_style smtbq

pair_style smtbq command
========================

Syntax
""""""


.. code-block:: LAMMPS

   pair_style smtbq

Examples
""""""""


.. code-block:: LAMMPS

   pair_style smtbq
   pair_coeff * * ffield.smtbq.Al2O3 O Al

Description
"""""""""""

This pair style computes a variable charge SMTB-Q (Second-Moment
tight-Binding QEq) potential as described in :ref:`SMTB-Q\_1 <SMTB-Q_1>` and
:ref:`SMTB-Q\_2 <SMTB-Q_2>`. Briefly, the energy of metallic-oxygen systems
is given by three contributions:

.. math::

   E_{tot} & =  E_{ES} + E_{OO} + E_{MO} \\
   E_{ES}  & =  \sum_i{\biggl[ \chi_{i}^{0}Q_i + \frac{1}{2}J_{i}^{0}Q_{i}^{2} +
   \frac{1}{2} \sum_{j\neq i}{ J_{ij}(r_{ij})f_{cut}^{R_{coul}}(r_{ij})Q_i Q_j } \biggr] } \\
   E_{OO}  & =  \sum_{i,j}^{i,j = O}{\biggl[Cexp( -\frac{r_{ij}}{\rho} ) - Df_{cut}^{r_1^{OO}r_2^{OO}}(r_{ij}) exp(Br_{ij})\biggr]}  \\
   E_{MO}  & =  \sum_i{E_{cov}^{i} + \sum_{j\neq i}{ Af_{cut}^{r_{c1}r_{c2}}(r_{ij})exp\bigl[-p(\frac{r_{ij}}{r_0} -1) \bigr] } }


where :math:`E_{tot}` is the total potential energy of the system,
:math:`E_{ES}` is the electrostatic part of the total energy,
:math:`E_{OO}` is the interaction between oxygen atoms and
:math:`E_{MO}` is a short-range interaction between metal and oxygen
atoms. This interactions depend on interatomic distance :math:`r_{ij}`
and/or the charge :math:`Q_{i}` of atoms *i*\ . Cut-off function enables
smooth convergence to zero interaction.

The parameters appearing in the upper expressions are set in the
ffield.SMTBQ.Syst file where Syst corresponds to the selected system
(e.g. field.SMTBQ.Al2O3). Examples for TiO2,
Al2O3 are provided.  A single pair\_coeff command
is used with the SMTBQ styles which provides the path to the potential
file with parameters for needed elements. These are mapped to LAMMPS
atom types by specifying additional arguments after the potential
filename in the pair\_coeff command. Note that atom type 1 must always
correspond to oxygen atoms. As an example, to simulate a TiO2 system,
atom type 1 has to be oxygen and atom type 2 Ti. The following
pair\_coeff command should then be used:


.. code-block:: LAMMPS

   pair_coeff * * PathToLammps/potentials/ffield.smtbq.TiO2 O Ti

The electrostatic part of the energy consists of two components 

self-energy of atom *i* in the form of a second order charge dependent
polynomial and a long-range Coulombic electrostatic interaction. The
latter uses the wolf summation method described in :ref:`Wolf <Wolf2>`,
spherically truncated at a longer cutoff, :math:`R_{coul}`. The
charge of each ion is modeled by an orbital Slater which depends on
the principal quantum number (\ *n*\ ) of the outer orbital shared by the
ion.

Interaction between oxygen, :math:`E_{OO}`, consists of two parts,
an attractive and a repulsive part. The attractive part is effective
only at short range (< :math:`r_2^{OO}`). The attractive
contribution was optimized to study surfaces reconstruction
(e.g. :ref:`SMTB-Q\_2 <SMTB-Q_2>` in TiO2) and is not necessary
for oxide bulk modeling. The repulsive part is the Pauli interaction
between the electron clouds of oxygen. The Pauli repulsion and the
coulombic electrostatic interaction have same cut off value. In the
ffield.SMTBQ.Syst, the keyword *'buck'* allows to consider only the
repulsive O-O interactions. The keyword *'buckPlusAttr'* allows to
consider the repulsive and the attractive O-O interactions.

The short-range interaction between metal-oxygen, :math:`E_{MO}` is
based on the second moment approximation of the density of states with
a N-body potential for the band energy term,
:math:`E^i_{cov}`, and a Born-Mayer type repulsive terms
as indicated by the keyword *'second\_moment'* in the
ffield.SMTBQ.Syst. The energy band term is given by:

.. math::

   E_{cov}^{i(i=M,O)} & = - \biggl\{\eta_i(\mu \xi^{0})^2 f_{cut}^{r_{c1}r_{c2}}(r_{ij})
   \biggl( \sum_{j(j=O,M)}{ exp[ -2q(\frac{r_{ij}}{r_0} - 1)] } \biggr) 
   \delta Q_i \bigl( 2\frac{n_0}{\eta_i} - \delta Q_i \bigr) \biggr\}^{1/2} \\
   \delta Q_i & =  | Q_i^{F} | - | Q_i |


where :math:\eta_i` is the stoichiometry of atom *i*\ ,
:math:`\delta Q_i` is the charge delocalization of atom *i*\ ,
compared to its formal charge
:math:`Q^F_i`. :math:`n_0`, the number of hybridized
orbitals, is calculated with to the atomic orbitals shared
:math:`d_i` and the stoichiometry
:math:`\eta_i`. :math:`r_{c1}` and :math:`r_{c2}` are the two
cutoff radius around the fourth neighbors in the cutoff function.

In the formalism used here, :math:`\xi^0` is the energy
parameter. :math:`\xi^0` is in tight-binding approximation the
hopping integral between the hybridized orbitals of the cation and the
anion. In the literature we find many ways to write the hopping
integral depending on whether one takes the point of view of the anion
or cation. These are equivalent vision. The correspondence between the
two visions is explained in appendix A of the article in the
SrTiO3 :ref:`SMTB-Q\_3 <SMTB-Q_3>` (parameter :math:`\beta` shown in
this article is in fact the :math:`\beta_O`). To summarize the
relationship between the hopping integral :math:`\xi^O`  and the
others, we have in an oxide CnOm the following
relationship:

.. math::

   \xi^0 & = \frac{\xi_O}{m} = \frac{\xi_C}{n} \\
   \frac{\beta_O}{\sqrt{m}} & = \frac{\beta_C}{\sqrt{n}} = \xi^0 \frac{\sqrt{m}+\sqrt{n}}{2}


Thus parameter :math:`\mu`, indicated above, is given by :math:`\mu = (\sqrt{n} + \sqrt{m}) / 2`

The potential offers the possibility to consider the polarizability of
the electron clouds of oxygen by changing the slater radius of the
charge density around the oxygen atoms through the parameters *rBB, rB and
rS* in the ffield.SMTBQ.Syst. This change in radius is performed
according to the method developed by E. Maras
:ref:`SMTB-Q\_2 <SMTB-Q_2>`. This method needs to determine the number of
nearest neighbors around the oxygen. This calculation is based on
first (:math:`r_{1n}`) and second (:math:`r_{2n}`) distances
neighbors.

The SMTB-Q potential is a variable charge potential. The equilibrium
charge on each atom is calculated by the electronegativity
equalization (QEq) method. See :ref:`Rick <Rick3>` for further detail. One
can adjust the frequency, the maximum number of iterative loop and the
convergence of the equilibrium charge calculation. To obtain the
energy conservation in NVE thermodynamic ensemble, we recommend to use
a convergence parameter in the interval 10e-5 -
10e-6 eV.

The ffield.SMTBQ.Syst files are provided for few systems. They consist
of nine parts and the lines beginning with '#' are comments (note that
the number of comment lines matter). The first sections are on the
potential parameters and others are on the simulation options and
might be modified. Keywords are character type and must be enclosed in
quotation marks ('').

1) Number of different element in the oxide:

* N_elem= 2 or 3
* Dividing line

2) Atomic parameters

For the anion (oxygen) 

* Name of element (char) and stoichiometry in oxide
* Formal charge and mass of element
* Principal quantum number of outer orbital n), electronegativity (:math:`\xi^0_i`) and hardness (:math:`J^0_i`)
* Ionic radius parameters  : max coordination number (\ *coordBB* = 6 by default), bulk coordination number *(coordB)*\ , surface coordination number  *(coordS)* and *rBB, rB and rS*  the slater radius for each coordination number. (**note : If you don't want to change the slater radius, use three identical radius values**)
* Number of orbital shared by the element in the oxide (:math:`d_i`)
* Dividing line

For each cations (metal):

* Name of element (char) and stoichiometry in oxide
* Formal charge and mass of element
* Number of electron in outer orbital *(ne)*\ , electronegativity (\ *&#967<sup>0</sup><sub>i</simulationub>*\ ), hardness (\ *J<sup>0</sup><sub>i</sub>*\ ) and *r<sub>Salter</sub>* the slater radius for the cation.
* Number of orbitals shared by the elements in the oxide (\ *d<sub>i</sub>*\ )
* Dividing line

3) Potential parameters:

* Keyword for element1, element2 and interaction potential ('second\_moment' or 'buck' or 'buckPlusAttr') between element 1 and 2.  If the potential is 'second\_moment', specify 'oxide' or 'metal' for metal-oxygen or metal-metal interactions respectively.
* Potential parameter: <pre><br/> If type of potential is 'second\_moment' : *A (eV)*\ , *p*\ , *&#958<sup>0</sup>* (eV) and *q* <br/> *r<sub>c1</sub>* (&#197), *r<sub>c2</sub>* (&#197) and *r<sub>0</sub>* (&#197) <br/> If type of potential is 'buck' : *C* (eV) and *&#961* (&#197) <br/> If type of potential is 'buckPlusAttr' : *C* (eV) and *&#961* (&#197) <br/> *D* (eV), *B* (&#197<sup>-1</sup>), *r<sub>1</sub><sup>OO</sup>* (&#197) and *r<sub>2</sub><sup>OO</sup>* (&#197) </pre>
* Dividing line

4) Tables parameters:

* Cutoff radius for the Coulomb interaction (\ *R<sub>coul</sub>*\ )
* Starting radius  (\ *r<sub>min</sub>* = 1,18845 &#197) and increments (\ *dr* = 0,001 &#197) for creating the potential table.
* Dividing line

5) Rick model parameter:

* *Nevery* : parameter to set the frequency (\ *1/Nevery*\ ) of the charge resolution. The charges are evaluated each *Nevery* time steps.
* Max number of iterative loop (\ *loopmax*\ ) and precision criterion (\ *prec*\ ) in eV of the charge resolution
* Dividing line

6) Coordination parameter:

* First (\ *r<sub>1n</sub>*\ ) and second (\ *r<sub>2n</sub>*\ ) neighbor distances in &#197
* Dividing line

7) Charge initialization mode:

* Keyword (\ *QInitMode*\ ) and initial oxygen charge (\ *Q<sub>init</sub>*\ ). If keyword = 'true', all oxygen charges are initially set equal to *Q<sub>init</sub>*\ . The charges on the cations are initially set in order to respect the neutrality of the box. If keyword = 'false', all atom charges are initially set equal to 0 if you use "create\_atom"#create\_atom command or the charge specified in the file structure using :doc:`read_data <read_data>` command.
* Dividing line

8) Mode for the electronegativity equalization (Qeq) 

* Keyword mode: <pre> <br/> QEqAll  (one QEq group) \|   no parameters <br/> QEqAllParallel (several QEq groups) \|   no parameters <br/> Surface \|   zlim   (QEq only for z>zlim)   </pre>
* Parameter if necessary
* Dividing line

9) Verbose 

* If you want the code to work in verbose mode or not : 'true' or 'false'
* If you want to print or not in file 'Energy\_component.txt' the three main contributions to the energy of the system according to the description presented above : 'true' or 'false' and *N<sub>Energy</sub>*\ . This option writes in file every *N<sub>Energy</sub>* time step. If the value is 'false' then *N<sub>Energy</sub>* = 0. The file take into account the possibility to have several QEq group *g* then it writes: time step, number of atoms in group *g*\ , electrostatic part of energy, *E<sub>ES</sub>*\ , the interaction between oxygen, *E<sub>OO</sub>*\ , and short range metal-oxygen interaction, *E<sub>MO</sub>*\ .
* If you want to print in file 'Electroneg\_component.txt' the electronegativity component (\ *&#8706E<sub>tot</sub> &#8260&#8706Q<sub>i</sub>*\ ) or not: 'true' or 'false' and *N<sub>Electroneg</sub>*\ .This option writes in file every *N<sub>Electroneg</sub>* time step. If the value is 'false' then *N<sub>Electroneg</sub>* = 0.  The file consist in atom number *i*\ , atom type (1 for oxygen and # higher than 1 for metal), atom position: *x*\ , *y* and *z*\ , atomic charge of atom *i*\ , electrostatic part of atom *i* electronegativity, covalent part of atom *i* electronegativity, the hopping integral of atom *i* *(Z&#946<sup>2</sup>)<sub>i<sub>* and box electronegativity.

.. note::

   This last option slows down the calculation dramatically.  Use
   only with a single processor simulation.


----------


**Mixing, shift, table, tail correction, restart, rRESPA info:**

This pair style does not support the :doc:`pair_modify <pair_modify>`
mix, shift, table, and tail options.

This pair style does not write its information to :doc:`binary restart files <restart>`, since it is stored in potential files.  Thus, you
needs to re-specify the pair\_style and pair\_coeff commands in an input
script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  It does not support the
*inner*\ , *middle*\ , *outer* keywords.


----------


**Restriction:**

This pair style is part of the USER-SMTBQ package and is only enabled
if LAMMPS is built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

This potential requires using atom type 1 for oxygen and atom type
higher than 1 for metal atoms.

This pair style requires the :doc:`newton <newton>` setting to be "on"
for pair interactions.

The SMTB-Q potential files provided with LAMMPS (see the potentials
directory) are parameterized for metal :doc:`units <units>`.


----------


**Citing this work:**

Please cite related publication: N. Salles, O. Politano, E. Amzallag
and R. Tetot, Comput. Mater. Sci. 111 (2016) 181-189


----------


.. _SMTB-Q\_1:



**(SMTB-Q\_1)** N. Salles, O. Politano, E. Amzallag, R. Tetot,
Comput. Mater. Sci. 111 (2016) 181-189

.. _SMTB-Q\_2:



**(SMTB-Q\_2)** E. Maras, N. Salles, R. Tetot, T. Ala-Nissila,
H. Jonsson, J. Phys. Chem. C 2015, 119, 10391-10399

.. _SMTB-Q\_3:



**(SMTB-Q\_3)** R. Tetot, N. Salles, S. Landron, E. Amzallag, Surface
Science 616, 19-8722 28 (2013)

.. _Wolf2:



**(Wolf)** D. Wolf, P. Keblinski, S. R. Phillpot, J. Eggebrecht, J Chem
Phys, 110, 8254 (1999).

.. _Rick3:



**(Rick)** S. W. Rick, S. J. Stuart, B. J. Berne, J Chem Phys 101, 6141
(1994).
