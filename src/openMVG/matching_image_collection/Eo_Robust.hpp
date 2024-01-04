// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Pierre MOULON, Romain JANVIER.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MATCHING_IMAGE_COLLECTION_EO_ROBUST_HPP
#define OPENMVG_MATCHING_IMAGE_COLLECTION_EO_ROBUST_HPP

#include <limits>
#include <utility>
#include <vector>

#include "openMVG/cameras/Camera_Pinhole.hpp"
#include "openMVG/matching/indMatch.hpp"
#include "openMVG/matching_image_collection/Geometric_Filter_utils.hpp"
#include "openMVG/multiview/solver_essential_kernel.hpp"
#include "openMVG/multiview/essential.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansac.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansacKernelAdaptator.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/types.hpp"

namespace openMVG {

namespace sfm {
struct Regions_Provider;
}

namespace matching_image_collection {

//-- A contrario essential matrix estimation template functor used for filter pair of putative correspondences
struct GeometricFilter_EOMatrix_RA
{
    GeometricFilter_EOMatrix_RA
    (
      double dPrecision = std::numeric_limits<double>::infinity(),
      uint32_t iteration = 1024
    ):
      m_dPrecision(dPrecision),
      m_stIteration(iteration),
      m_E(Mat3::Identity())
    {
    }

    /// Robust fitting according an Orthographic Epipolar Geometry criteria
    template<typename Regions_or_Features_ProviderT>
    bool Robust_estimation
    (
      const sfm::SfM_Data * sfm_data,
      const std::shared_ptr<Regions_or_Features_ProviderT> & regions_provider,
      const Pair pairIndex,
      const matching::IndMatches & vec_PutativeMatches,
      matching::IndMatches & geometric_inliers
    )
    {
      geometric_inliers.clear();

      // Get back corresponding view index
      const IndexT
        iIndex = pairIndex.first,
        jIndex = pairIndex.second;

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
      if (!isPinhole(cam_I->getType()) || !isPinhole(cam_J->getType()))
        return false;

      //--
      // Get corresponding point regions arrays
      //--

      Mat2X xI,xJ;
      MatchesPairToMat(pairIndex, vec_PutativeMatches, sfm_data, regions_provider, xI, xJ);

      // Update precision if required (normalization from image to camera plane):
      if (m_dPrecision != std::numeric_limits<double>::infinity())
      {
        m_dPrecision = (cam_I->imagePlane_toCameraPlaneError(Square(m_dPrecision)) +
                        cam_J->imagePlane_toCameraPlaneError(Square(m_dPrecision))) / 2.;
      }

      //--
      // Robust estimation
      //--

      // Define the Kernel
      // --- using KernelType = essential::kernel::ThreePointKernel;
      using KernelType =
        robust::ACKernelAdaptorEssentialOrtho<
          essential::kernel::ThreePointKernel,
          essential::kernel::OrthographicSymmetricEpipolarDistanceError,
          Mat3>;

      const KernelType kernel(
        (*cam_I)(xI),
        sfm_data->GetViews().at(iIndex)->ui_width,
        sfm_data->GetViews().at(iIndex)->ui_height,
        (*cam_J)(xJ),
        sfm_data->GetViews().at(jIndex)->ui_width,
        sfm_data->GetViews().at(jIndex)->ui_height
      );

      // Robustly estimate the model with AC-RANSAC
      std::vector<uint32_t> vec_inliers;

      const auto ACRansacOut = ACRANSAC(
        kernel, vec_inliers, m_stIteration, &m_E, m_dPrecision);

      if (vec_inliers.size() > KernelType::MINIMUM_SAMPLES * 2.5)
      {
        // update geometric_inliers
        geometric_inliers.reserve(vec_inliers.size());
        for (const uint32_t & index : vec_inliers)
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

  uint32_t m_stIteration; // maximal number of iteration for robust estimation
  double m_dPrecision;    // upper_bound precision used for robust estimation
  //
  //-- Stored data
  Mat3 m_E;
};

} //namespace matching_image_collection
} // namespace openMVG

#endif //OPENMVG_MATCHING_IMAGE_COLLECTION_EO_ROBUST_HPP
