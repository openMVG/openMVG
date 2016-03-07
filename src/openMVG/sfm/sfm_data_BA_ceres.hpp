// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_DATA_BA_CERES_HPP
#define OPENMVG_SFM_DATA_BA_CERES_HPP

#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_BA.hpp"
#include "openMVG/sfm/sfm_data_BA_ceres_camera_functor.hpp"
#include "ceres/ceres.h"

namespace openMVG {
namespace sfm {

/// Create the appropriate cost functor according the provided input camera intrinsic model
ceres::CostFunction * IntrinsicsToCostFunction(
  cameras::IntrinsicBase * intrinsic,
  const Vec2 & observation);

class Bundle_Adjustment_Ceres : public Bundle_Adjustment
{
  public:
  struct BA_options
  {
    bool bVerbose_;
    unsigned int nb_threads_;
    bool bCeres_summary_;
    ceres::LinearSolverType linear_solver_type_;
    ceres::PreconditionerType preconditioner_type_;
    ceres::SparseLinearAlgebraLibraryType sparse_linear_algebra_library_type_;

    BA_options(const bool bVerbose = true, bool bmultithreaded = true);
  };
  private:
    BA_options openMVG_options_;

  public:
  Bundle_Adjustment_Ceres(Bundle_Adjustment_Ceres::BA_options options = BA_options());

  bool Adjust
  (
    // the SfM scene to refine
    SfM_Data & sfm_data,
    // tell which parameter needs to be adjusted
    Parameter_Adjustment_Option adjustment_option = Parameter_Adjustment_Option::ADJUST_ALL
  );  // tell if the structure will be refined
};

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_DATA_BA_CERES_HPP
