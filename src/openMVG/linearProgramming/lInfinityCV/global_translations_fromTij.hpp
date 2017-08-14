// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_LINFINITY_COMPUTER_VISION_GLOBAL_TRANSLATIONS_FROMBEARING_HPP
#define OPENMVG_LINFINITY_COMPUTER_VISION_GLOBAL_TRANSLATIONS_FROMBEARING_HPP

#include <utility>
#include <vector>

#include "openMVG/linearProgramming/linearProgrammingInterface.hpp"
#include "openMVG/multiview/translation_averaging_common.hpp"

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
void EncodeTi_from_tij
(
  const size_t nTranslation,
  const std::vector<openMVG::RelativeInfo_Vec> & vec_relative_motion_groups,
  sRMat & A, Vec & C,
  std::vector<LP_Constraints::eLP_SIGN> & vec_sign,
  std::vector<double> & vec_costs,
  std::vector< std::pair<double,double> > & vec_bounds
);

//-- Estimate the translation from heading relative translations of triplets.
//- Translation directions must not be normalized (in this way relative scale
//-  of relative motions is kept and colinear motion is supported).
struct Tifromtij_ConstraintBuilder
{
  explicit Tifromtij_ConstraintBuilder
  (
    const std::vector< openMVG::RelativeInfo_Vec > & vec_relative
  );

  /// Setup constraints for the global translations problem,
  ///  in the LP_Constraints_Sparse object.
  bool Build(LP_Constraints_Sparse & constraint);

  // Internal data
  size_t Ncam_;
  const std::vector< openMVG::RelativeInfo_Vec > & vec_relative_; // /!\ memory Alias
};

} // namespace lInfinityCV
} // namespace openMVG

#endif // OPENMVG_LINFINITY_COMPUTER_VISION_GLOBAL_TRANSLATIONS_FROMBEARING_HPP
