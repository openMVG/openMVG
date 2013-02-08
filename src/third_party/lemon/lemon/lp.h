/* -*- mode: C++; indent-tabs-mode: nil; -*-
 *
 * This file is a part of LEMON, a generic C++ optimization library.
 *
 * Copyright (C) 2003-2010
 * Egervary Jeno Kombinatorikus Optimalizalasi Kutatocsoport
 * (Egervary Research Group on Combinatorial Optimization, EGRES).
 *
 * Permission to use, modify and distribute this software is granted
 * provided that this copyright notice appears in all copies. For
 * precise terms see the accompanying LICENSE file.
 *
 * This software is provided "AS IS" with no warranty of any kind,
 * express or implied, and with no claim as to its suitability for any
 * purpose.
 *
 */

#ifndef LEMON_LP_H
#define LEMON_LP_H

#include<lemon/config.h>


#ifdef LEMON_HAVE_GLPK
#include <lemon/glpk.h>
#elif LEMON_HAVE_CPLEX
#include <lemon/cplex.h>
#elif LEMON_HAVE_SOPLEX
#include <lemon/soplex.h>
#elif LEMON_HAVE_CLP
#include <lemon/clp.h>
#endif

///\file
///\brief Defines a default LP solver
///\ingroup lp_group
namespace lemon {

#ifdef DOXYGEN
  ///The default LP solver identifier

  ///The default LP solver identifier.
  ///\ingroup lp_group
  ///
  ///Currently, the possible values are \c GLPK, \c CPLEX,
  ///\c SOPLEX or \c CLP
#define LEMON_DEFAULT_LP SOLVER
  ///The default LP solver

  ///The default LP solver.
  ///\ingroup lp_group
  ///
  ///Currently, it is either \c GlpkLp, \c CplexLp, \c SoplexLp or \c ClpLp
  typedef GlpkLp Lp;

  ///The default MIP solver identifier

  ///The default MIP solver identifier.
  ///\ingroup lp_group
  ///
  ///Currently, the possible values are \c GLPK or \c CPLEX
#define LEMON_DEFAULT_MIP SOLVER
  ///The default MIP solver.

  ///The default MIP solver.
  ///\ingroup lp_group
  ///
  ///Currently, it is either \c GlpkMip or \c CplexMip
  typedef GlpkMip Mip;
#else
#ifdef LEMON_HAVE_GLPK
# define LEMON_DEFAULT_LP GLPK
  typedef GlpkLp Lp;
# define LEMON_DEFAULT_MIP GLPK
  typedef GlpkMip Mip;
#elif LEMON_HAVE_CPLEX
# define LEMON_DEFAULT_LP CPLEX
  typedef CplexLp Lp;
# define LEMON_DEFAULT_MIP CPLEX
  typedef CplexMip Mip;
#elif LEMON_HAVE_SOPLEX
# define DEFAULT_LP SOPLEX
  typedef SoplexLp Lp;
#elif LEMON_HAVE_CLP
# define DEFAULT_LP CLP
  typedef ClpLp Lp;
#endif
#endif

} //namespace lemon

#endif //LEMON_LP_H
