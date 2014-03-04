########################################################
L infinity solvers for computer vision
########################################################

openMVG propose Linear programming based solver for various problem in computer vision by minimizing of the maximal error the residuals between observations and estimated parameters (The L infinity norm).

openMVG implements problems introduced by [LinfNorm]_ and generalized by [LinfNormGeneric].

Rather than considering quadratic constraints that require SOCP (Second Orde Cone Programming) we consider their LP (linear program) equivalent. It makes usage of residual error expressed with absolute error ( ``|a|<b`` inequality is transformed in two linear inequalities ``a<b`` and ``-b<-a``. It makes the solving faster and constraint easier to express (see. [Arnak]_ for more explanation).

- N-view triangulation [LinfNorm]_,
- Resection or pose matrix estimation [LinfNorm]_,
- Estimation of translations and structure from known rotations,

  - two formulation are implemented,

    - the simple one [LinfNorm]_,
    - the robust based on slack variables [OlssonDuality]_.

- Registration of relative translations to compute global translations [ICCV13]_,

  - using triplets of translations.



.. [LinfNorm] L infinity Minimization in Geometric Reconstruction Problems. Richard I. Hartley, Frederik Schaffalitzky. CVPR 2004.

.. [LinfNormGeneric] "Multiple-View Geometry under the L infty Norm." Authors: Fredrik Kahl, Richard Hartley. 2008. and "Multiple View Geometry and the L infty -norm". Fredrik Kahl. ICCV 2005.

.. [Arnak] Robust estimation for an inverse problem arising in multiview geometry. Arnak S. Dalalyan, Renaud Keriven. In J. Math. Imaging Vision, 2012.


.. [OlssonDuality] Outlier Removal Using Duality. Carl Olsson, Anders Eriksson and Richard Hartley, Richard. CVPR 2010.


.. [ICCV13] Global Fusion of Relative Motions for Robust, Accurate and Scalable Structure from Motion. Pierre Moulon, Pascal Monasse and Renaud Marlet. ICCV 2013.
