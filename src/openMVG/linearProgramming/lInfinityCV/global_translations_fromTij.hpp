// Copyright (c) 2012 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_LINFINITY_COMPUTER_VISION_GLOBAL_TRANSLATIONS_FROMBEARING_H_
#define OPENMVG_LINFINITY_COMPUTER_VISION_GLOBAL_TRANSLATIONS_FROMBEARING_H_

#include "openMVG/numeric/numeric.h"
#include "openMVG/multiview/translation_averaging_common.hpp"
#include "openMVG/linearProgramming/linearProgrammingInterface.hpp"
#include <fstream>
#include <utility>
#include <vector>
#include <set>

#ifdef _MSC_VER
#pragma warning( once : 4267 ) //warning C4267: 'argument' : conversion from 'size_t' to 'const int', possible loss of data
#endif

//------------------
//-- Bibliography --
//------------------
//- [1] "Global Fusion of Relative Motions for Robust, Accurate and Scalable Structure from Motion."
//- Authors: Pierre MOULON, Pascal MONASSE and Renaud MARLET.
//- Date: December 2013.
//- Conference: ICCV.

namespace openMVG   {
namespace lInfinityCV  {

using namespace linearProgramming;

// Setup the linear program to solve the union of trifocal tensors heading
//  directions in a common global coordinate system.
// This implementation is a generalization of the LINEAR PROGRAM (9) page 5 of [1]
// -> This implementation can deal with groups of relative motions.
//    You can mix bearing vectors of 2-view, 3-view, X-view configuration.
//    Each group will have it's own shared scaling factor.
//--
inline void EncodeTi_from_tij(
    const size_t nTranslation,
    const std::vector<openMVG::RelativeInfo_Vec> & vec_relative_motion_groups,
    sRMat & A, Vec & C,
    std::vector<LP_Constraints::eLP_SIGN> & vec_sign,
    std::vector<double> & vec_costs,
    std::vector< std::pair<double,double> > & vec_bounds)
{
  // Build Constraint matrix.

  const size_t Ncam = (size_t) nTranslation;

  size_t Nrelative = 0;
  for (const openMVG::RelativeInfo_Vec & iter : vec_relative_motion_groups)
  {
    Nrelative += iter.size();
  }

  const size_t transStart  = 0;
  const size_t lambdaStart = 3 * Ncam;
  const size_t gammaStart = lambdaStart + vec_relative_motion_groups.size();

#undef TVAR
#undef LAMBDAVAR
#undef GAMMAVAR
# define TVAR(i, el) (transStart + 3*(i) + (el))  // translation (X,Y,Z)
# define LAMBDAVAR(j) (lambdaStart + (int)(j)) // One per relative group of translations
# define GAMMAVAR gammaStart

  const size_t Nconstraint = Nrelative * 6;
  const size_t NVar =
    // GlobalTranslation(XYZ) + OneScaleFactorPerMotionGroup + Gamma
    3 * Ncam + vec_relative_motion_groups.size() + 1;

  A.resize(Nconstraint, NVar);

  C.resize(Nconstraint, 1);
  C.fill(0.0);
  vec_sign.resize(Nconstraint);

  // By default set free variables:
  vec_bounds = std::vector< std::pair<double,double> >(NVar);
  fill( vec_bounds.begin(), vec_bounds.end(),
    std::make_pair((double)-1e+30, (double)1e+30));
  //  Make first camera at origin (translation ambiguity)
  vec_bounds[TVAR(0,0)].first = vec_bounds[TVAR(0,0)].second = 0;
  vec_bounds[TVAR(0,1)].first = vec_bounds[TVAR(0,1)].second = 0;
  vec_bounds[TVAR(0,2)].first = vec_bounds[TVAR(0,2)].second = 0;
  // Make lambda variables between 1 and large number => constraint that lambda_ij > 1
  for (size_t k = 0; k < vec_relative_motion_groups.size(); ++k)
    vec_bounds[lambdaStart + k].first = 1;

  // Setup gamma >= 0
  vec_bounds[vec_bounds.size()-1].first = 0.0;

  //-- Minimize gamma
  vec_costs.resize(NVar);
  std::fill(vec_costs.begin(), vec_costs.end(), 0.0);
  vec_costs[GAMMAVAR] = 1.0;
  //--

  size_t rowPos = 0;
  size_t group_id = 0;
  for (const openMVG::RelativeInfo_Vec & iter : vec_relative_motion_groups)
  {
    for (const relativeInfo & rel : iter)
    {
      const size_t i   = rel.first.first;
      const size_t j   = rel.first.second;
      const Mat3 & Rij = rel.second.first;
      const Vec3 & tij = rel.second.second;

      // | T_j - R_ij T_i - Lambda_ij t_ij | < Gamma
      // Absolute constraint transformed in two sign constraints
      //   T_j - R_ij T_i - Lambda_ij t_ij < Gamma
      //   T_j - R_ij T_i - Lambda_ij t_ij > - Gamma

      // For X, Y, Z axis:
      for (int l = 0; l < 3; ++l)
      {
        // T_j
        A.coeffRef(rowPos, TVAR(j, l)) = 1;

        //- R_ij T_i
        A.coeffRef(rowPos, TVAR(i, 0)) = - Rij(l, 0);
        A.coeffRef(rowPos, TVAR(i, 1)) = - Rij(l, 1);
        A.coeffRef(rowPos, TVAR(i, 2)) = - Rij(l, 2);

        // - Lambda_ij t_ij
        A.coeffRef(rowPos, LAMBDAVAR(group_id)) = - tij(l);

        // - gamma
        A.coeffRef(rowPos, GAMMAVAR) = -1;

        // < 0
        vec_sign[rowPos] = LP_Constraints::LP_LESS_OR_EQUAL;
        C(rowPos) = 0;
        ++rowPos;

        // ---------
        // Opposite constraint
        //   T_j - R_ij T_i - Lambda_ij t_ij > - Gamma
        // ---------

        // T_j
        A.coeffRef(rowPos, TVAR(j, l)) = 1;

        //- R_ij T_i
        A.coeffRef(rowPos, TVAR(i, 0)) = - Rij(l, 0);
        A.coeffRef(rowPos, TVAR(i, 1)) = - Rij(l, 1);
        A.coeffRef(rowPos, TVAR(i, 2)) = - Rij(l, 2);

        // - Lambda_ij t_ij
        A.coeffRef(rowPos, LAMBDAVAR(group_id)) = - tij(l);

        // + gamma
        A.coeffRef(rowPos, GAMMAVAR) = 1;

        // > 0
        vec_sign[rowPos] = LP_Constraints::LP_GREATER_OR_EQUAL;
        C(rowPos) = 0;
        ++rowPos;
      }
    }
    ++group_id;
  }
#undef TVAR
#undef LAMBDAVAR
#undef GAMMAVAR
}

//-- Estimate the translation from heading relative translations of triplets.
//- Translation directions must not be normalized (in this way relative scale
//-  of relative motions is kept and colinear motion is supported).
struct Tifromtij_ConstraintBuilder
{
  Tifromtij_ConstraintBuilder(
    const std::vector< openMVG::RelativeInfo_Vec > & vec_relative
  )
  : vec_relative_(vec_relative)
  {
    //Count the number of camera that are represented
    std::set<size_t> countSet;
    for (const openMVG::RelativeInfo_Vec & iter : vec_relative_)
    {
      for (const relativeInfo & rel : iter)
      {
        countSet.insert(rel.first.first);
        countSet.insert(rel.first.second);
      }
    }
    Ncam_ = countSet.size();
  }

  /// Setup constraints for the global translations problem,
  ///  in the LP_Constraints_Sparse object.
  bool Build(LP_Constraints_Sparse & constraint)
  {
    EncodeTi_from_tij(
      Ncam_,
      vec_relative_,
      constraint.constraint_mat_,
      constraint.constraint_objective_,
      constraint.vec_sign_,
      constraint.vec_cost_,
      constraint.vec_bounds_);

    // it's a minimization problem over the gamma variable
    constraint.bminimize_ = true;

    //-- Setup additional information about the Linear Program constraint.
    // We look for :
    //  - #translations parameters,
    //  - #relative lambda factors (one per group of relative motions),
    //  - one gamma parameter.
    constraint.nbParams_ = Ncam_ * 3 + vec_relative_.size() + 1;
    return true;
  }

  // Internal data
  size_t Ncam_;
  const std::vector< openMVG::RelativeInfo_Vec > & vec_relative_; // /!\ memory Alias
};

} // namespace lInfinityCV
} // namespace openMVG

#endif // OPENMVG_LINFINITY_COMPUTER_VISION_GLOBAL_TRANSLATIONS_FROMBEARING_H_

