/* -*- mode: C++; indent-tabs-mode: nil; -*-
 *
 * This file is a part of LEMON, a generic C++ optimization library.
 *
 * Copyright (C) 2003-2013
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


#if LEMON_DEFAULT_LP == LEMON_GLPK_ || LEMON_DEFAULT_MIP == LEMON_GLPK_
#include <lemon/glpk.h>
#endif
#if LEMON_DEFAULT_LP == LEMON_CPLEX_ || LEMON_DEFAULT_MIP == LEMON_CPLEX_
#include <lemon/cplex.h>
#endif
#if LEMON_DEFAULT_LP == LEMON_SOPLEX_
#include <lemon/soplex.h>
#endif
#if LEMON_DEFAULT_LP == LEMON_CLP_
#include <lemon/clp.h>
#endif
#if LEMON_DEFAULT_MIP == LEMON_CBC_
#include <lemon/cbc.h>
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
  ///Currently, the possible values are \c LEMON_GLPK_, \c LEMON_CPLEX_,
  ///\c LEMON_SOPLEX_ or \c LEMON_CLP_
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
  ///Currently, the possible values are \c LEMON_GLPK_, \c LEMON_CPLEX_
  ///or \c LEMON_CBC_
#define LEMON_DEFAULT_MIP SOLVER
  ///The default MIP solver.

  ///The default MIP solver.
  ///\ingroup lp_group
  ///
  ///Currently, it is either \c GlpkMip, \c CplexMip , \c CbcMip
  typedef GlpkMip Mip;
#else
#if LEMON_DEFAULT_LP == LEMON_GLPK_
  typedef GlpkLp Lp;
#elif LEMON_DEFAULT_LP == LEMON_CPLEX_
  typedef CplexLp Lp;
#elif LEMON_DEFAULT_LP == LEMON_SOPLEX_
  typedef SoplexLp Lp;
#elif LEMON_DEFAULT_LP == LEMON_CLP_
  typedef ClpLp Lp;
#endif
#if LEMON_DEFAULT_MIP == LEMON_GLPK_
  typedef GlpkMip Mip;
#elif LEMON_DEFAULT_MIP == LEMON_CPLEX_
  typedef CplexMip Mip;
#elif LEMON_DEFAULT_MIP == LEMON_CBC_
  typedef CbcMip Mip;
#endif
#endif

} //namespace lemon

#endif //LEMON_LP_H
