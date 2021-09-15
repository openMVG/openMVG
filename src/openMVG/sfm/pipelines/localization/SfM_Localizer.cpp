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
#include "openMVG/multiview/solver_resection_up2p_kukelova.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_BA.hpp"
#include "openMVG/sfm/sfm_data_BA_ceres.hpp"
#include "openMVG/sfm/sfm_landmark.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansac.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansacKernelAdaptator.hpp"

#include <memory>
#include <utility>

namespace openMVG
{
/// Pose/Resection Kernel adapter for the A contrario model estimator with
///  known camera intrinsics.
template <typename SolverArg,
  typename ModelArg = Mat34>
class ACKernelAdaptorResection_Intrinsics
{
public:
  using Solver = SolverArg;
  using Model = ModelArg;

  ACKernelAdaptorResection_Intrinsics
  (
    const Mat & x2d, // Undistorted 2d feature_point location
    const Mat & x3D, // 3D corresponding points
    const cameras::IntrinsicBase * camera
  ):x2d_(x2d),
    x3D_(x3D),
    logalpha0_(log10(M_PI)),
    N1_(Mat3::Identity()),
    camera_(camera)
  {
    N1_.diagonal().head(2) *= camera->imagePlane_toCameraPlaneError(1.0);
    assert(2 == x2d_.rows());
    assert(3 == x3D_.rows());
    assert(x2d_.cols() == x3D_.cols());
    bearing_vectors_= camera->operator()(x2d_);
  }

  enum { MINIMUM_SAMPLES = Solver::MINIMUM_SAMPLES };
  enum { MAX_MODELS = Solver::MAX_MODELS };

  void Fit(const std::vector<uint32_t> &samples, std::vector<Model> *models) const {
    Solver::Solve(ExtractColumns(bearing_vectors_, samples), // bearing vectors
                  ExtractColumns(x3D_, samples), // 3D points
                  models); // Found model hypothesis
  }

  void Errors(const Model & model, std::vector<double> & vec_errors) const
  {
    // Convert the found model into a Pose3
    const Vec3 t = model.block(0, 3, 3, 1);
    const geometry::Pose3 pose(model.block(0, 0, 3, 3),
                               - model.block(0, 0, 3, 3).transpose() * t);

    vec_errors.resize(x2d_.cols());

    const bool ignore_distortion = true; // We ignore distortion since we are using undistorted bearing vector as input

    for (Mat::Index sample = 0; sample < x2d_.cols(); ++sample)
    {
      vec_errors[sample] = (camera_->residual(pose(x3D_.col(sample)),
                              x2d_.col(sample),
                              ignore_distortion) * N1_(0,0)).squaredNorm();
    }
  }

  size_t NumSamples() const { return x2d_.cols(); }

  void Unnormalize(Model * model) const {
  }

  double logalpha0() const {return logalpha0_;}
  double multError() const {return 1.0;} // point to point error
  Mat3 normalizer1() const {return Mat3::Identity();}
  Mat3 normalizer2() const {return N1_;}
  double unormalizeError(double val) const {return sqrt(val) / N1_(0,0);}

private:
  Mat x2d_, bearing_vectors_;
  const Mat & x3D_;
  Mat3 N1_;
  double logalpha0_;  // Alpha0 is used to make the error adaptive to the image size
  const cameras::IntrinsicBase * camera_;   // Intrinsic camera parameter
};
} // namespace openMVG

namespace openMVG {
namespace sfm {

  bool SfM_Localizer::Localize
  (
    const resection::SolverType & solver_type,
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

    switch (solver_type)
    {
      case resection::SolverType::DLT_6POINTS:
      {
        //--
        // Classic resection (try to compute the entire P matrix)
        using SolverType = openMVG::resection::kernel::SixPointResectionSolver;
        MINIMUM_SAMPLES = SolverType::MINIMUM_SAMPLES;

        using KernelType =
          openMVG::robust::ACKernelAdaptorResection<
            SolverType,
            resection::SquaredPixelReprojectionError,
            openMVG::robust::UnnormalizerResection,
            Mat34>;

        KernelType kernel(resection_data.pt2D, image_size.first, image_size.second,
          resection_data.pt3D);
        // Robust estimation of the pose and its precision
        const std::pair<double,double> ACRansacOut =
          openMVG::robust::ACRANSAC(kernel,
                                    resection_data.vec_inliers,
                                    resection_data.max_iteration,
                                    &P,
                                    dPrecision,
                                    true);
        // Update the upper bound precision of the model found by AC-RANSAC
        resection_data.error_max = ACRansacOut.first;
      }
      break;
      case resection::SolverType::P3P_NORDBERG_ECCV18:
      {
        if (!optional_intrinsics)
        {
          std::cerr << "Intrinsic data is required for P3P solvers." << std::endl;
          return false;
        }
        //--
        // Since the intrinsic data is known, compute only the pose
        using SolverType = openMVG::euclidean_resection::P3PSolver_Nordberg;
        MINIMUM_SAMPLES = SolverType::MINIMUM_SAMPLES;

        using KernelType =
          ACKernelAdaptorResection_Intrinsics<
            SolverType,
            Mat34>;

        KernelType kernel(resection_data.pt2D, resection_data.pt3D, optional_intrinsics);
        // Robust estimation of the pose matrix and its precision
        const auto ACRansacOut =
          openMVG::robust::ACRANSAC(kernel,
                                    resection_data.vec_inliers,
                                    resection_data.max_iteration,
                                    &P,
                                    dPrecision,
                                    true);
        // Update the upper bound precision of the model found by AC-RANSAC
        resection_data.error_max = ACRansacOut.first;
      }
      break;
      case resection::SolverType::P3P_KE_CVPR17:
      {
        if (!optional_intrinsics)
        {
          std::cerr << "Intrinsic data is required for P3P solvers." << std::endl;
          return false;
        }
        //--
        // Since the intrinsic data is known, compute only the pose
        using SolverType = openMVG::euclidean_resection::P3PSolver_Ke;
        MINIMUM_SAMPLES = SolverType::MINIMUM_SAMPLES;

        using KernelType =
          ACKernelAdaptorResection_Intrinsics<
            SolverType,
            Mat34>;

        KernelType kernel(resection_data.pt2D, resection_data.pt3D, optional_intrinsics);
        // Robust estimation of the pose matrix and its precision
        const auto ACRansacOut =
          openMVG::robust::ACRANSAC(kernel,
                                    resection_data.vec_inliers,
                                    resection_data.max_iteration,
                                    &P,
                                    dPrecision,
                                    true);
        // Update the upper bound precision of the model found by AC-RANSAC
        resection_data.error_max = ACRansacOut.first;
      }
      break;
      case resection::SolverType::P3P_KNEIP_CVPR11:
      {
        if (!optional_intrinsics)
        {
          std::cerr << "Intrinsic data is required for P3P solvers." << std::endl;
          return false;
        }
        //--
        // Since the intrinsic data is known, compute only the pose
        using SolverType = openMVG::euclidean_resection::P3PSolver_Kneip;
        MINIMUM_SAMPLES = SolverType::MINIMUM_SAMPLES;

        using KernelType =
          ACKernelAdaptorResection_Intrinsics<
            SolverType,
            Mat34>;

        KernelType kernel(resection_data.pt2D, resection_data.pt3D, optional_intrinsics);

        // Robust estimation of the pose matrix and its precision
        const auto ACRansacOut =
          openMVG::robust::ACRANSAC(kernel,
                                    resection_data.vec_inliers,
                                    resection_data.max_iteration,
                                    &P,
                                    dPrecision,
                                    true);
        // Update the upper bound precision of the model found by AC-RANSAC
        resection_data.error_max = ACRansacOut.first;
      }
      break;
      case resection::SolverType::UP2P_KUKELOVA_ACCV10:
      {
        if (!optional_intrinsics)
        {
          std::cerr << "Intrinsic data is required for P3P solvers." << std::endl;
          return false;
        }
        //--
        // Since the intrinsic data is known, compute only the pose
        using SolverType = openMVG::euclidean_resection::UP2PSolver_Kukelova;
        MINIMUM_SAMPLES = SolverType::MINIMUM_SAMPLES;

        using KernelType =
          ACKernelAdaptorResection_Intrinsics<
            SolverType,
            Mat34>;

        KernelType kernel(resection_data.pt2D, resection_data.pt3D, optional_intrinsics);

        // Robust estimation of the pose matrix and its precision
        const auto ACRansacOut =
          openMVG::robust::ACRANSAC(kernel,
                                    resection_data.vec_inliers,
                                    resection_data.max_iteration,
                                    &P,
                                    dPrecision,
                                    true);
        // Update the upper bound precision of the model found by AC-RANSAC
        resection_data.error_max = ACRansacOut.first;
      }
      break;
      default:
      {
        std::cerr << "Unknown absolute pose solver type." << std::endl;
        return false;
      }
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
        intrinsics->updateFromParams(shared_intrinsics->getParams());
    }

    return b_BA_Status;
  }

} // namespace sfm
} // namespace openMVG
