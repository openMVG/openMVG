============
Ceres Solver
============

Ceres Solver [#f1]_ is an open source C++ library for modeling and
solving large, complicated optimization problems. It can be used to
solve `Non-linear Least Squares`_ problems with bounds constraints and
general unconstrained optimization problems. It is a mature, feature
rich, and performant library that has been used in production at
Google since 2010. For more, see :doc:`features`.

`ceres-solver@googlegroups.com
<https://groups.google.com/forum/?fromgroups#!forum/ceres-solver>`_ is
the place for discussions and questions about Ceres Solver. We use the
`GitHub Issue Tracker
<https://github.com/ceres-solver/ceres-solver/issues>`_ to manage bug
reports and feature requests.


.. toctree::
   :maxdepth: 1
   :hidden:

   features
   installation
   tutorial
   derivatives
   nnls_modeling
   nnls_solving
   nnls_covariance
   gradient_solver
   faqs
   users
   contributing
   version_history
   bibliography
   license

.. _Non-linear Least Squares: http://en.wikipedia.org/wiki/Non-linear_least_squares


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
