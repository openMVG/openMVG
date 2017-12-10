// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MATCHING_IMAGE_COLLECTION_E_SPHERICAL_ACROBUST_ANGULAR_HPP
#define OPENMVG_MATCHING_IMAGE_COLLECTION_E_SPHERICAL_ACROBUST_ANGULAR_HPP

#include <limits>
#include <utility>
#include <vector>

#include "openMVG/cameras/Camera_Intrinsics.hpp"
#include "openMVG/matching/indMatch.hpp"
#include "openMVG/matching_image_collection/Geometric_Filter_utils.hpp"
#include "openMVG/multiview/essential.hpp"
#include "openMVG/multiview/motion_from_essential.hpp"
#include "openMVG/multiview/solver_essential_eight_point.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansac.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansacKernelAdaptator.hpp"
#include "openMVG/sfm/pipelines/sfm_robust_model_estimation.hpp"
#include "openMVG/sfm/sfm_data.hpp"

namespace openMVG { namespace sfm { struct Regions_Provider; } }

namespace openMVG {
namespace matching_image_collection {

//-- A contrario essential matrix estimation template functor used for filter pair of putative correspondences
struct GeometricFilter_ESphericalMatrix_AC_Angular
{
  GeometricFilter_ESphericalMatrix_AC_Angular(
    double precision_upper_bound = std::numeric_limits<double>::infinity(),
    size_t iteration = 1024)
    : m_precision_upper_bound(precision_upper_bound),
      m_stIteration(iteration),
      m_E(Mat3::Identity()),
      m_precision_upper_bound_robust(std::numeric_limits<double>::infinity())
  {

  };

  /// Robust fitting of the ESSENTIAL matrix
  template<typename Regions_or_Features_ProviderT>
  bool Robust_estimation(
    const sfm::SfM_Data * sfm_data,
    const std::shared_ptr<Regions_or_Features_ProviderT> & regions_provider,
    const Pair pairIndex,
    const matching::IndMatches & vec_PutativeMatches,
    matching::IndMatches & geometric_inliers)
  {
    using namespace openMVG;
    using namespace openMVG::robust;
    geometric_inliers.clear();

    // Get back corresponding view index
    const IndexT iIndex = pairIndex.first;
    const IndexT jIndex = pairIndex.second;

    //--
    // Reject pair with missing Intrinsic information
    //--

    const sfm::View
      * view_I = sfm_data->views.at(iIndex).get(),
      * view_J = sfm_data->views.at(jIndex).get();

     // Check that valid cameras can be retrieved for the pair of views
    const cameras::IntrinsicBase
      * cam_I =
        sfm_data->GetIntrinsics().count(view_I->id_intrinsic) ?
          sfm_data->GetIntrinsics().at(view_I->id_intrinsic).get() : nullptr,
      * cam_J =
        sfm_data->GetIntrinsics().count(view_J->id_intrinsic) ?
          sfm_data->GetIntrinsics().at(view_J->id_intrinsic).get() : nullptr;

    if (!cam_I || !cam_J)
      return false;

    //--
    // Get corresponding point regions arrays
    //--

    Mat2X xI,xJ;
    MatchesPairToMat(pairIndex, vec_PutativeMatches, sfm_data, regions_provider, xI, xJ);

    const Mat3X
      xI_bearing_vector = (*cam_I)(xI),
      xJ_bearing_vector = (*cam_J)(xJ);

    //--
    // Robust estimation
    //--

    // Define the AContrario angular error adaptor
    typedef openMVG::robust::ACKernelAdaptor_AngularRadianError<
        // Use the 8 point solver in order to estimate E
        openMVG::EightPointRelativePoseSolver,
        openMVG::AngularError,
        Mat3>
        KernelType;

    KernelType kernel(xI_bearing_vector, xJ_bearing_vector);

    // Robustly estimate the Essential matrix with A Contrario ransac
    const double upper_bound_precision =
     (m_precision_upper_bound != std::numeric_limits<double>::infinity())?
        D2R(m_precision_upper_bound) : std::numeric_limits<double>::infinity();
    std::vector<uint32_t> vec_inliers;
    const auto ac_ransac_output =
      ACRANSAC(kernel, vec_inliers, m_stIteration, &m_E, upper_bound_precision);

    const double & threshold = ac_ransac_output.first;

    // Decompose the essential matrix and keep the best solution (if any)
    {
      geometry::Pose3 relative_pose;
      std::vector<uint32_t> inliers_indexes;
      std::vector<Vec3> inliers_X;
      if (RelativePoseFromEssential(
        xI_bearing_vector, xJ_bearing_vector,
        m_E, vec_inliers,
        &relative_pose, &inliers_indexes, &inliers_X))
      {
        m_relativePose = relative_pose;
        vec_inliers = inliers_indexes;
      }
      else
      {
        vec_inliers.clear();
      }
    }
    if (vec_inliers.size() > KernelType::MINIMUM_SAMPLES * 2.5)
    {
      m_precision_upper_bound_robust = ac_ransac_output.first;
      // update geometric_inliers
      geometric_inliers.reserve(vec_inliers.size());
      for (const size_t & index : vec_inliers)
      {
        geometric_inliers.push_back( vec_PutativeMatches[index] );
      }
      return true;
    }
    else
    {
      vec_inliers.clear();
      return false;
    }
  }

  bool Geometry_guided_matching
  (
    const sfm::SfM_Data * sfm_data,
    const std::shared_ptr<sfm::Regions_Provider> & regions_provider,
    const Pair pairIndex,
    const double dDistanceRatio,
    matching::IndMatches & matches
  )
  {
    return false;
  }

  double m_precision_upper_bound;  // upper_bound precision used for robust estimation
  size_t m_stIteration; // maximal number of iteration for robust estimation
  //
  //-- Stored data
  Mat3 m_E;
  geometry::Pose3 m_relativePose;
  double m_precision_upper_bound_robust;
};

} // namespace matching_image_collection
} // namespace openMVG

#endif // OPENMVG_MATCHING_IMAGE_COLLECTION_E_SPHERICAL_ACROBUST_ANGULAR_HPP
