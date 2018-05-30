// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/multiview/translation_averaging_common.hpp"
#include "openMVG/multiview/translation_averaging_solver.hpp"
#include "openMVG/numeric/eigen_alias_definition.hpp"
#include "openMVG/types.hpp"

#ifdef OPENMVG_USE_OPENMP
#include <omp.h>
#endif

#include <ceres/ceres.h>
#include <ceres/rotation.h>

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
  explicit SmallScaleError
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
  const std::vector<openMVG::RelativeInfo_Vec > & vec_relative_group_estimates,
  std::vector<Eigen::Vector3d> & translations,
  const double d_l1_loss_threshold
)
{
  //-- Count:
  //- #poses are used by the relative position estimates
  //- #relative estimates we will use

  std::set<unsigned int> count_set;
  unsigned int relative_info_count = 0;
  for (const openMVG::RelativeInfo_Vec & iter : vec_relative_group_estimates)
  {
    for (const relativeInfo & it_relative_motion : iter)
    {
      ++relative_info_count;
      count_set.insert(it_relative_motion.first.first);
      count_set.insert(it_relative_motion.first.second);
    }
  }
  const IndexT nb_poses = count_set.size();

  //--
  // Build the parameters arrays:
  //--
  // - camera translations
  // - relative translation scales (one per group)
  // - relative rotations

  std::vector<double> vec_translations(relative_info_count*3, 1.0);
  const unsigned nb_scales = vec_relative_group_estimates.size();
  std::vector<double> vec_scales(nb_scales, 1.0);

  // Setup the relative rotations array (angle axis parametrization)
  std::vector<double> vec_relative_rotations(relative_info_count*3, 0.0);
  unsigned int cpt = 0;
  for (const openMVG::RelativeInfo_Vec & iter : vec_relative_group_estimates)
  {
    for (const relativeInfo & info : iter)
    {
      ceres::RotationMatrixToAngleAxis(
        (const double*)info.second.first.data(),
        &vec_relative_rotations[cpt]);
      cpt += 3;
    }
  }

  ceres::LossFunction * loss =
    (d_l1_loss_threshold < 0) ? nullptr : new ceres::SoftLOneLoss(d_l1_loss_threshold);

  // Add constraints to the minimization problem
  ceres::Problem problem;
  //
  // A. Add cost functors:
  cpt = 0;
  IndexT scale_idx = 0;
  for (const openMVG::RelativeInfo_Vec & iter : vec_relative_group_estimates)
  {
    for (const relativeInfo & info : iter)
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
        &vec_relative_rotations[cpt],
        &vec_scales[scale_idx]);
      // the relative rotation is set as constant
      problem.SetParameterBlockConstant(&vec_relative_rotations[cpt]);
      cpt+=3;
    }
    ++scale_idx; // One scale per relative_motion group
  }

  // B. Add constraint over the scale factors:
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
  if (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::SUITE_SPARSE))
  {
    options.sparse_linear_algebra_library_type = ceres::SUITE_SPARSE;
    options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
  }
  else if (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::CX_SPARSE))
  {
    options.sparse_linear_algebra_library_type = ceres::CX_SPARSE;
    options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
  }
  else if (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::EIGEN_SPARSE))
  {
    options.sparse_linear_algebra_library_type = ceres::EIGEN_SPARSE;
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
#if CERES_VERSION_MAJOR < 2
  options.num_linear_solver_threads = omp_get_max_threads();
#endif
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
  for (unsigned int i = 0; i < nb_poses; ++i, cpt+=3)
  {
    translations[i] << vec_translations[cpt], vec_translations[cpt+1], vec_translations[cpt+2];
  }
  return true;
}

} // namespace openMVG
