README
======

:Authors: Juan Luis Cano Rodríguez <jl.cano@alumnos.upm.es>,
    José Luis Cercós Pita <jl.cercos@upm.es>

Overview
--------

.. image:: http://unmaintained.tech/badge.svg
   :target: http://unmaintained.tech/
   :alt: No Maintenance Intended 

This is the code for the final master project "Integration of a free
finite element method structural solver with an SPH code for the study
of fluid-structure interaction".

* ``lib``: Library code containing reusable classes and functions

  - ``test``: Unit tests corresponding to the library code

* ``examples``: Sample programs to demonstrate how to use the library code
* ``report``: Code corresponding to the figures and results of the report
* ``buildscripts``: Build scripts (conda recipes) to easily install the
  dependencies

.. warning:: This code is unfinished and the project has been abandoned. Use
  at your own risk.

Requirements
------------

The code requires FEniCS 1.6.0 to run properly. Notice that it is not
backwards compatible with previous versions. An ``environment.yml`` file is
provided::

  $ conda env create

Demos
-----

To run the demos, use::

  $ ./run.sh --run

