// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_DATA_BA_CERES_HPP
#define OPENMVG_SFM_DATA_BA_CERES_HPP

#include "openMVG/sfm/sfm_data_BA.hpp"
#include "openMVG/numeric/numeric.h"
#include "ceres/types.h"
#include "ceres/cost_function.h"

namespace openMVG {

namespace cameras{
class IntrinsicBase;
}

namespace sfm {

struct SfM_Data;

/// Create the appropriate cost functor according the provided input camera intrinsic model
/// Can be residual cost functor can be weighetd if desired (default 0.0 means no weight).
ceres::CostFunction * IntrinsicsToCostFunction
(
  cameras::IntrinsicBase * intrinsic,
  const Vec2 & observation,
  const double weight = 0.0
);

class Bundle_Adjustment_Ceres : public Bundle_Adjustment
{
  public:
  struct BA_Ceres_options
  {
    bool bVerbose_;
    unsigned int nb_threads_;
    bool bCeres_summary_;
    ceres::LinearSolverType linear_solver_type_;
    ceres::PreconditionerType preconditioner_type_;
    ceres::SparseLinearAlgebraLibraryType sparse_linear_algebra_library_type_;
    double parameter_tolerance_;
    bool bUse_loss_function_;

    BA_Ceres_options(const bool bVerbose = true, bool bmultithreaded = true);
  };
  private:
    BA_Ceres_options ceres_options_;

  public:
  Bundle_Adjustment_Ceres(Bundle_Adjustment_Ceres::BA_Ceres_options options = BA_Ceres_options());

  BA_Ceres_options & ceres_options();

  bool Adjust
  (
    // the SfM scene to refine
    SfM_Data & sfm_data,
    // tell which parameter needs to be adjusted
    const Optimize_Options options
  ) override;
};

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_DATA_BA_CERES_HPP
