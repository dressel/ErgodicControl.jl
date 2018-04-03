.. ErgodicControl.jl documentation master file, created by
   sphinx-quickstart on Wed Feb  8 16:36:49 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

ErgodicControl.jl - Design of Ergodic Trajectories
=====================================================

ErgodicControl.jl is a Julia package for performing ergodic control---that is, designing trajectories that are ergodic with respect to some spatial distribution.
Currently, two trajectory generation methods are implemented: projection-based optimization and spectral multi-scale coverage.
The plan is to implement others, including sequential action control, HEDAC, and sampling based methods.

As of April 2018, this code is under active development and the interface may change.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

    Installation <installation.rst>
    Ergodic Control Overview <ergodic_control.rst>
    Basic Types <basic_types.rst>
    Ergodic Manager <ergodic_manager.rst>
    Trajectory Manager <trajectory_manager.rst>
    Generating Trajectories <generating_trajectories.rst>
    Visuals <visuals.rst>
    Examples <examples.rst>



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
