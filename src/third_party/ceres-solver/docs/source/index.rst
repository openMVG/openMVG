.. Ceres Solver documentation master file, created by
   sphinx-quickstart on Sat Jan 19 00:07:33 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

============
Ceres Solver
============

.. toctree::
   :maxdepth: 3
   :hidden:

   features
   building
   tutorial
   modeling
   solving
   faqs
   contributing
   version_history
   history
   bibliography
   license

Ceres Solver is an open source C++ library for modeling and solving
large complicated `nonlinear least squares`_ problems. It is a feature
rich, mature and performant library which has been used in production
since 2010. At Google, Ceres Solver is used to:

* Estimate the pose of `Street View`_ cars, aircrafts, and satellites.
* Build 3D models for `PhotoTours`_.
* Estmate satellite image sensor characteristics.
* Stitch `panoramas`_ or apply `Lens Blur`_ on Android.
* Solve `bundle adjustment`_ and SLAM problems in `Project Tango`_.

Outside Google, Ceres is used for solving problems in computer vision,
computer graphics, astronomy and physics. e.g., `Willow Garage`_ uses
it to solve SLAM problems and `Blender`_ uses it for for planar
tracking and bundle adjustment.

.. _nonlinear least squares: http://en.wikipedia.org/wiki/Non-linear_least_squares
.. _fitting curves: http://en.wikipedia.org/wiki/Nonlinear_regression
.. _bundle adjustment: http://en.wikipedia.org/wiki/Structure_from_motion
.. _Street View: http://youtu.be/z00ORu4bU-A
.. _PhotoTours: http://google-latlong.blogspot.com/2012/04/visit-global-landmarks-with-photo-tours.html
.. _panoramas: http://www.google.com/maps/about/contribute/photosphere/
.. _Project Tango: https://www.google.com/atap/projecttango/
.. _Blender: http://mango.blender.org/development/planar-tracking-preview/
.. _Willow Garage: https://www.willowgarage.com/blog/2013/08/09/enabling-robots-see-better-through-improved-camera-calibration
.. _Lens Blur: http://googleresearch.blogspot.com/2014/04/lens-blur-in-new-google-camera-app.html

Getting started
---------------

* Download the `latest stable release
  <http://ceres-solver.org/ceres-solver-1.9.0.tar.gz>`_ or clone the
  Git repository for the latest development version.

  .. code-block:: bash

       git clone https://ceres-solver.googlesource.com/ceres-solver

* Read the :ref:`chapter-tutorial`, browse the chapters on the
  :ref:`chapter-modeling` API and the :ref:`chapter-solving` API.
* Join the `mailing list
  <https://groups.google.com/forum/?fromgroups#!forum/ceres-solver>`_
  and ask questions.
* File bugs, feature requests in the `issue tracker
  <https://code.google.com/p/ceres-solver/issues/list>`_.


Cite Us
-------
If you use Ceres Solver for a publication, please cite it as::

    @misc{ceres-solver,
      author = "Sameer Agarwal and Keir Mierle and Others",
      title = "Ceres Solver",
      howpublished = "\url{http://ceres-solver.org}",
    }
