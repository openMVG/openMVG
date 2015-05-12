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
   api
   faqs
   users
   contributing
   version_history
   bibliography
   license

Ceres Solver [#f1]_ is an open source C++ library for modeling and
solving large, complicated optimization problems. It is a feature
rich, mature and performant library which has been used in production
at Google since 2010. Ceres Solver can solve two kinds of problems.

1. `Non-linear Least Squares`_ problems with bounds constraints.
2. General unconstrained optimization problems.

.. _Non-linear Least Squares: http://en.wikipedia.org/wiki/Non-linear_least_squares

Getting started
===============

* Download the `latest stable release
  <http://ceres-solver.org/ceres-solver-1.10.0.tar.gz>`_ or clone the
  Git repository for the latest development version.

  .. code-block:: bash

       git clone https://ceres-solver.googlesource.com/ceres-solver

* Read the :ref:`chapter-tutorial` and browse the :ref:`chapter-api`.
* Join the `mailing list
  <https://groups.google.com/forum/?fromgroups#!forum/ceres-solver>`_
  and ask questions.
* File bugs, feature requests in the `issue tracker
  <https://code.google.com/p/ceres-solver/issues/list>`_.


Cite Us
=======

If you use Ceres Solver for a publication, please cite it as::

    @misc{ceres-solver,
      author = "Sameer Agarwal and Keir Mierle and Others",
      title = "Ceres Solver",
      howpublished = "\url{http://ceres-solver.org}",
    }


.. rubric:: Footnotes

.. [#f1] While there is some debate as to who invented the method of
         Least Squares [Stigler]_, there is no questioning the fact
         that it was `Carl Friedrich Gauss
         <http://www-groups.dcs.st-and.ac.uk/~history/Biographies/Gauss.html>`_
         who brought it to the attention of the world. Using just 22
         observations of the newly discovered asteroid `Ceres
         <http://en.wikipedia.org/wiki/Ceres_(dwarf_planet)>`_, Gauss
         used the method of least squares to correctly predict when
         and where the asteroid will emerge from behind the Sun
         [TenenbaumDirector]_. We named our solver after Ceres to
         celebrate this seminal event in the history of astronomy,
         statistics and optimization.
