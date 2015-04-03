TODO
----

.. important::
    Put the right amount of effort into this project. Not a bit less,
    not a bit more.

Items in this TODO:

a) Things about the document or the report itself: things to explain or clarify -> These are low priority
b) Things about the code organization, cleaning of the setup -> These are top priority
c) Things about expanding the work and moving forward -> These are medium priority

.. note:: Some things need an Internet connection

Top priority
~~~~~~~~~~~~

Work in progress
~~~~~~~~~~~~~~~~

* Where do I put the client-server architecture?
* Do any patches to AQUAgpusph remain? Probably the VTK linking problem
  still stands, you will need inspiration with DOLFIN code itself to produce
  a proper CMake file

Other
~~~~~

* Redo the aquaclient.py inspired by FEniCS.py
* Redo the figures with the TiKZ backend
* Sad status of CBC.Solve: flow is unmaintained, rock is missing, beat is
  undocumented, twist mixes tabs and spaces, common does not even pass
  static analysis
* Make up your mind about Neumann and Dirichlet conditions
* Chapter 29 can be of general interest
* The icing in the cake: explain the presence of possible hydroelastic
  effects and applications of this work, given that in the dam-break
  case these won't be noticed
* Include report
* Produce proper tests and measure coverage, I will need the py_func attribute
* Make a different repo for both? Probably at a later stage

  a) Library + Examples + Unit tests
  b) Report with the LaTeX source too

* Client/Server: A REQ/REP architecture has been chosen and it seems good
  enough for our purposes with the process we designed
* Write a chapter on how did you do the demo to clarify the architecture of AQUAgpusph
* Yak shaving: don't forget to stop your DigitalOcean machine if you are
  not using it

Process
~~~~~~~

* run.sh clears replaced files, calls create_sensors.py and then Create.py
  and starts the AQUAgpusph binary with Main.xml as input
* create_sensors.py reads the FEniCS mesh and creates Sensors.dat with
  one point for each DOF
* Create.py sets the problem parameters, creates the geometry, takes the
  template files and replace the values

