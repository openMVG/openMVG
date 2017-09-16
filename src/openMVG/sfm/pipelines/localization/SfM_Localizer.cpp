// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/pipelines/localization/SfM_Localizer.hpp"

#include "openMVG/cameras/Camera_Common.hpp"
#include "openMVG/cameras/Camera_Intrinsics.hpp"
#include "openMVG/cameras/Camera_Pinhole.hpp"
#include "openMVG/multiview/solver_resection_kernel.hpp"
#include "openMVG/multiview/solver_resection_p3p.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_BA.hpp"
#include "openMVG/sfm/sfm_data_BA_ceres.hpp"
#include "openMVG/sfm/sfm_landmark.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansac.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansacKernelAdaptator.hpp"

#include <memory>
#include <utility>

namespace openMVG {
namespace sfm {

  struct ResectionSquaredResidualError {
    // Compute the residual of the projection distance(pt2D, Project(P,pt3D))
    // Return the squared error
    static double Error(const Mat34 & P, const Vec2 & pt2D, const Vec3 & pt3D){
      const Vec2 x = Project(P, pt3D);
      return (x - pt2D).squaredNorm();
    }
  };

  bool SfM_Localizer::Localize
  (
    const Pair & image_size,
    const cameras::IntrinsicBase * optional_intrinsics,
    Image_Localizer_Match_Data & resection_data,
    geometry::Pose3 & pose
  )
  {
    // --
    // Compute the camera pose (resectioning)
    // --
    Mat34 P;
    resection_data.vec_inliers.clear();

    // Setup the admissible upper bound residual error
    const double dPrecision =
      resection_data.error_max == std::numeric_limits<double>::infinity() ?
      std::numeric_limits<double>::infinity() :
      Square(resection_data.error_max);

    size_t MINIMUM_SAMPLES = 0;
    const cameras::Pinhole_Intrinsic * pinhole_cam = dynamic_cast<const cameras::Pinhole_Intrinsic *>(optional_intrinsics);
    if (pinhole_cam == nullptr)
    {
      //--
      // Classic resection (try to compute the entire P matrix)
      using SolverType = openMVG::resection::kernel::SixPointResectionSolver;
      MINIMUM_SAMPLES = SolverType::MINIMUM_SAMPLES;

      using KernelType =
        openMVG::robust::ACKernelAdaptorResection<
        SolverType,
        ResectionSquaredResidualError,
        openMVG::robust::UnnormalizerResection,
        Mat34>;

      KernelType kernel(resection_data.pt2D, image_size.first, image_size.second,
        resection_data.pt3D);
      // Robust estimation of the Projection matrix and it's precision
      const std::pair<double,double> ACRansacOut =
        openMVG::robust::ACRANSAC(kernel, resection_data.vec_inliers, resection_data.max_iteration, &P, dPrecision, true);
      // Update the upper bound precision of the model found by AC-RANSAC
      resection_data.error_max = ACRansacOut.first;
    }
    else
    {
      //--
      // Since K calibration matrix is known, compute only [R|t]
      using SolverType = openMVG::euclidean_resection::P3PSolver;
      MINIMUM_SAMPLES = SolverType::MINIMUM_SAMPLES;

      using KernelType =
        openMVG::robust::ACKernelAdaptorResection_K<
          SolverType,
          ResectionSquaredResidualError,
          Mat34>;

      KernelType kernel(resection_data.pt2D, resection_data.pt3D, pinhole_cam->K());
      // Robust estimation of the Projection matrix and it's precision
      const std::pair<double,double> ACRansacOut =
        openMVG::robust::ACRANSAC(kernel, resection_data.vec_inliers, resection_data.max_iteration, &P, dPrecision, true);
      // Update the upper bound precision of the model found by AC-RANSAC
      resection_data.error_max = ACRansacOut.first;
    }

    // Test if the mode support some points (more than those required for estimation)
    const bool bResection = (resection_data.vec_inliers.size() > 2.5 * MINIMUM_SAMPLES);

    if (bResection)
    {
      resection_data.projection_matrix = P;
      Mat3 K, R;
      Vec3 t;
      KRt_From_P(P, &K, &R, &t);
      pose = geometry::Pose3(R, -R.transpose() * t);
    }

    std::cout << "\n"
      << "-------------------------------" << "\n"
      << "-- Robust Resection " << "\n"
      << "-- Resection status: " << bResection << "\n"
      << "-- #Points used for Resection: " << resection_data.pt2D.cols() << "\n"
      << "-- #Points validated by robust Resection: " << resection_data.vec_inliers.size() << "\n"
      << "-- Threshold: " << resection_data.error_max << "\n"
      << "-------------------------------" << std::endl;

    return bResection;
  }

  bool SfM_Localizer::RefinePose
  (
    cameras::IntrinsicBase * intrinsics,
    geometry::Pose3 & pose,
    Image_Localizer_Match_Data & matching_data,
    bool b_refine_pose,
    bool b_refine_intrinsic
  )
  {
    if (!b_refine_pose && !b_refine_intrinsic)
    {
      // Nothing to do (There is no parameter to refine)
      return false;
    }

    // Setup a tiny SfM scene with the corresponding 2D-3D data
    SfM_Data sfm_data;
    // view
    sfm_data.views.insert({0, std::make_shared<View>("",0, 0, 0)});
    // pose
    sfm_data.poses[0] = pose;
    // intrinsic
    std::shared_ptr<cameras::IntrinsicBase> shared_intrinsics(intrinsics->clone());
    sfm_data.intrinsics[0] = shared_intrinsics;
    // structure data (2D-3D correspondences)
    for (size_t i = 0; i < matching_data.vec_inliers.size(); ++i)
    {
      const size_t idx = matching_data.vec_inliers[i];
      Landmark landmark;
      landmark.X = matching_data.pt3D.col(idx);
      landmark.obs[0] = Observation(matching_data.pt2D.col(idx), UndefinedIndexT);
      sfm_data.structure[i] = std::move(landmark);
    }

    // Configure BA options (refine the intrinsic and the pose parameter only if requested)
    const Optimize_Options ba_refine_options
    (
      (b_refine_intrinsic) ? cameras::Intrinsic_Parameter_Type::ADJUST_ALL : cameras::Intrinsic_Parameter_Type::NONE,
      (b_refine_pose) ? Extrinsic_Parameter_Type::ADJUST_ALL : Extrinsic_Parameter_Type::NONE,
      Structure_Parameter_Type::NONE // STRUCTURE must remain constant
    );
    Bundle_Adjustment_Ceres bundle_adjustment_obj;
    const bool b_BA_Status = bundle_adjustment_obj.Adjust(
      sfm_data,
      ba_refine_options);
    if (b_BA_Status)
    {
      pose = sfm_data.poses[0];
      if (b_refine_intrinsic)
        *intrinsics = *shared_intrinsics;
    }

    return b_BA_Status;
  }

} // namespace sfm
} // namespace openMVG
