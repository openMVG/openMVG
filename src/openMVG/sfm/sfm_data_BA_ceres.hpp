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
    bool m_bVerbose;
    unsigned int m_nbThreads;
    bool _bCeres_Summary;
    ceres::LinearSolverType m_linear_solver_type;
    ceres::PreconditionerType m_preconditioner_type;
    ceres::SparseLinearAlgebraLibraryType m_sparse_linear_algebra_library_type;
    double m_parameter_tolerance;

    BA_options(const bool bVerbose = true, bool bmultithreaded = true);
  };
  private:
    BA_options m_openMVG_options;

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
