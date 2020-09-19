LAMMPS version |version| Documentation
######################################

LAMMPS stands for **L**\ arge-scale **A**\ tomic/**M**\ olecular
**M**\ assively **P**\ arallel **S**\ imulator.

LAMMPS is a classical molecular dynamics simulation code with a focus
on materials modeling.  It was designed to run efficiently on parallel
computers.  It was developed originally at Sandia National
Laboratories, a US Department of Energy facility.  The majority of
funding for LAMMPS has come from the US Department of Energy (DOE).
LAMMPS is an open-source code, distributed freely under the terms of
the GNU Public License (GPL).

The `LAMMPS website <lws_>`_ has a variety of information about the code.
It includes links to an on-line version of this manual, a `mailing list <https://lammps.sandia.gov/mail.html>`_ where users can post
questions, and a `GitHub site <https://github.com/lammps/lammps>`_ where
all LAMMPS development is coordinated.

----------

The content for this manual is part of the LAMMPS distribution.  You
can build a local copy of the Manual as HTML pages or a PDF file, by
following the steps on the :doc:`Manual build <Manual_build>` doc page.
The manual is organized in two parts:
1) A :ref:`User documentation <user_documentation>` for how to install
and use LAMMPS and 2) a :ref:`Programmer documentation <programmer_documentation>`
for how to write programs using the LAMMPS library or how to modify LAMMPS.

----------

Once you are familiar with LAMMPS, you may want to bookmark :doc:`this page <Commands_all>` since it gives quick access to a doc page for
every LAMMPS command.

.. _lws: https://lammps.sandia.gov

User Documentation
******************

..   :caption: User Documentation
.. _user_documentation:
.. toctree::
   :maxdepth: 2
   :numbered: 3
   :name: userdoc
   :includehidden:

   Intro
   Install
   Build
   Run_head
   Commands
   Packages
   Speed
   Howto
   Examples
   Tools
   Python_head
   Errors
   Manual_build

Programmer Documentation
************************
..   :caption: Programmer Documentation
.. _programmer_documentation:
.. toctree::
   :maxdepth: 2
   :numbered: 3
   :name: progdoc
   :includehidden:

   pg_library
   Modify
   pg_developer

.. toctree::
   :caption: Index
   :name: index
   :hidden:

   commands_list
   fixes
   computes
   pairs
   bonds
   angles
   dihedrals
   impropers
   fix_modify_atc_commands

Indices and tables
******************

.. only:: html

   * :ref:`genindex`
   * :ref:`search`
