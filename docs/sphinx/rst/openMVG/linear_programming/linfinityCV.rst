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

- Registration of relative translations to compute global translations [GlobalACSfM]_,

  - using triplets of translations.

