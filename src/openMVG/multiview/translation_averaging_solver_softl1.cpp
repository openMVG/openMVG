
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/multiview/translation_averaging_common.hpp"
#include "openMVG/multiview/translation_averaging_solver.hpp"
#include "openMVG/numeric/numeric.h"
#include "openMVG/types.hpp"

#include "ceres/ceres.h"
#include "ceres/rotation.h"

#include <vector>

namespace openMVG {

// Main Cost functor for translation averaging:
// measure the consistency (residual error) from the relative translations (constant) to the scales and global camera translations
struct RelativeTranslationError
{
  RelativeTranslationError
  (
    double t_ij_x,
    double t_ij_y,
    double t_ij_z
  )
  : t_ij_x_(t_ij_x), t_ij_y_(t_ij_y), t_ij_z_(t_ij_z)
  {}

  template <typename T>
  bool operator()
  (
    const T* const t_i,
    const T* const t_j,
    const T* const R_ij,
    const T* const s_ij,
    T* residuals
  ) const
  {
    //
    // residual = t_j - R_ij t_i - s_ij * t_ij
    //
    T rotated_t_i[3];
    // Rotate the point according the relative local rotation
    ceres::AngleAxisRotatePoint(R_ij, t_i, rotated_t_i);

    residuals[0] = T(t_j[0] - rotated_t_i[0] - *s_ij * t_ij_x_);
    residuals[1] = T(t_j[1] - rotated_t_i[1] - *s_ij * t_ij_y_);
    residuals[2] = T(t_j[2] - rotated_t_i[2] - *s_ij * t_ij_z_);

    return true;
  }

  double t_ij_x_, t_ij_y_, t_ij_z_;
};

// Cost penalizing scales smaller than 1.
struct SmallScaleError
{
  SmallScaleError
  (
    double weight = 1.0
  )
      : weight_(weight)
  {}

  template <typename T>
  bool operator()
  (
    const T* const s_ij, T* residual
  ) const
  {
    residual[0] = (*s_ij > T(1.0)) ? T(0.0) : (T(weight_) * (T(1.0) - *s_ij));
    return true;
  }

  double weight_;
};


bool solve_translations_problem_softl1
(
  const std::vector<openMVG::relativeInfo > & vec_initial_estimates,
  const bool b_translation_triplets,
  const int nb_poses,
  std::vector<Eigen::Vector3d> & translations,
  const double d_l1_loss_threshold
)
{
  ceres::Problem problem;

  // Build the parameters arrays:
  // - camera translation
  // - relative translation group scales
  // - relative rotations

  std::vector<double> vec_translations(3*nb_poses, 1.0);
  const unsigned nb_scales = vec_initial_estimates.size() / (b_translation_triplets ? 3 : 1);
  std::vector<double> vec_scales(nb_scales, 1.0);

  if (!b_translation_triplets)
  {
    // use random initialization, since using only single bearing vector results
    //  in a is less conditionned system.
    for (int i=0; i<nb_scales; ++i) {
      vec_scales[i] = (double)rand() / RAND_MAX;
    }

    for (int i=0; i<3*nb_poses; ++i) {
      vec_translations[i] = (double)rand() / RAND_MAX;
    }
  }

  // Relative rotations array
  std::vector<double> vec_relative_rotations(vec_initial_estimates.size()*3, 0.0);
  size_t cpt = 0;
  for (const openMVG::relativeInfo & info : vec_initial_estimates)
  {
    ceres::RotationMatrixToAngleAxis(
      (const double*)info.second.first.data(),
      &vec_relative_rotations[cpt]);
    cpt += 3;
  }

  ceres::LossFunction * loss =
    (d_l1_loss_threshold < 0) ? nullptr : new ceres::SoftLOneLoss(d_l1_loss_threshold);

  // Add constraints to the minimization
  //
  // A. Add cost functor from camera translation to the relative informations
  cpt = 0;
  IndexT scale_idx = 0;
  for (const openMVG::relativeInfo & info : vec_initial_estimates)
  {
    const Pair & ids = info.first;
    const IndexT I = ids.first;
    const IndexT J = ids.second;
    const Vec3 t_ij = info.second.second;

    // Each Residual block takes 2 camera translations & the relative rotation & a scale
    // and outputs a 3 dimensional residual.
    ceres::CostFunction* cost_function =
        new ceres::AutoDiffCostFunction<RelativeTranslationError, 3, 3, 3, 3,1>(
            new RelativeTranslationError(t_ij(0), t_ij(1), t_ij(2)));
    problem.AddResidualBlock(
      cost_function,
       loss,
       &vec_translations[I*3],
       &vec_translations[J*3],
       &vec_relative_rotations[cpt*3],
       &vec_scales[scale_idx]);
    // the relative rotation is set as constant
    problem.SetParameterBlockConstant(&vec_relative_rotations[cpt*3]);
    if (cpt % (b_translation_triplets ? 3 : 1) == 0 && cpt != 0)
      scale_idx += 1;
    ++cpt;
  }

  // B. Constraint the scale factors:
  //  Prefer scale > 1, since a trivial solution is translations = {0,...,0}).
  for (unsigned i = 0; i < nb_scales; ++i)
  {
    ceres::CostFunction* cost_function =
        new ceres::AutoDiffCostFunction<SmallScaleError, 1, 1>(
            new SmallScaleError(1.0));

    problem.AddResidualBlock(cost_function, nullptr, &vec_scales[i]);
  }
  // Set one center as known (to fix the gauge freedom)
  vec_translations[0] = vec_translations[1] = vec_translations[2] = 0.0;
  problem.SetParameterBlockConstant(&vec_translations[0]);

  // Solve
  ceres::Solver::Options options;
  options.minimizer_progress_to_stdout = false;
  if (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::SUITE_SPARSE) ||
      ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::CX_SPARSE) ||
      ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::EIGEN_SPARSE))
  {
    options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
  }
  else
  {
    options.linear_solver_type = ceres::DENSE_NORMAL_CHOLESKY;
  }
  options.max_num_iterations = std::max(50, (int)(nb_scales * 2));
  options.minimizer_progress_to_stdout = false;
  options.logging_type = ceres::SILENT;
#ifdef OPENMVG_USE_OPENMP
  options.num_threads = omp_get_max_threads();
  options.num_linear_solver_threads = omp_get_max_threads();
#endif // OPENMVG_USE_OPENMP
  ceres::Solver::Summary summary;
  ceres::Solve(options, &problem, &summary);

  if (!summary.IsSolutionUsable())
  {
    std::cout << summary.FullReport() << std::endl;
    return false;
  }

  // Fill the global translations array
  translations.resize(nb_poses);
  cpt = 0;
  for (unsigned i = 0; i < nb_poses; ++i, cpt+=3)
  {
    translations[i] << vec_translations[cpt], vec_translations[cpt+1], vec_translations[cpt+2];
  }
  return true;
}

} // namespace openMVG
