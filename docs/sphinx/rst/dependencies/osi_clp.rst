############################
OSI CLP
############################

========
USAGE
========

openMVG uses the [OSI]_ and the [CLP]_ solver in order to solve linear programs [LP]_.

CLP has been choosen because it is known to support problems of up to 1.5 million constraints [CLP FAQ]_.

[LPSOLVE]_ has been tested but tests shown that it is less reliable (sometimes, there is no convergence to a solution).

===========
Description
===========

[OSI]_
  Osi (Open Solver Interface) provides an abstract base class to a generic linear programming (LP) solver, along with derived classes for specific solvers. Many applications may be able to use the Osi to insulate themselves from a specific LP solver. (CLP and MOSEK backend are supported in openMVG).

[CLP]_
  Clp (Coin-or linear programming) is an open-source linear programming solver written in C++.

.. [OSI] https://projects.coin-or.org/Osi
.. [CLP] https://projects.coin-or.org/Clp/wiki
.. [CLP_FAQ] https://projects.coin-or.org/Clp/wiki/FAQ
.. [LPSOLVE] http://lpsolve.sourceforge.net/5.5/
