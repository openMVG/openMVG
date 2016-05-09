
// Copyright (c) 2014, 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "openMVG/multiview/solver_homography_kernel.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansac.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansacKernelAdaptator.hpp"
#include "openMVG/robust_estimation/guided_matching.hpp"

#include "openMVG/matching/indMatch.hpp"
#include "openMVG/matching/indMatchDecoratorXY.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/pipelines/sfm_regions_provider.hpp"
#include "openMVG/matching_image_collection/Geometric_Filter_utils.hpp"

namespace openMVG {
namespace matching_image_collection {

//-- A contrario homography matrix estimation template functor used for filter pair of putative correspondences
struct GeometricFilter_HMatrix_AC
{
  GeometricFilter_HMatrix_AC(
    double dPrecision = std::numeric_limits<double>::infinity(),
    size_t iteration = 1024)
    : m_dPrecision(dPrecision), m_stIteration(iteration), m_H(Mat3::Identity()),
      m_dPrecision_robust(std::numeric_limits<double>::infinity()){};

  /// Robust fitting of the HOMOGRAPHY matrix
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
    // Get corresponding point regions arrays
    //--

    Mat xI,xJ;
    MatchesPairToMat(pairIndex, vec_PutativeMatches, sfm_data, regions_provider, xI, xJ);

    //--
    // Robust estimation
    //--

    // Define the AContrario adapted Homography matrix solver
    typedef ACKernelAdaptor<
      openMVG::homography::kernel::FourPointSolver,
      openMVG::homography::kernel::AsymmetricError,
      UnnormalizerI,
      Mat3>
      KernelType;

    KernelType kernel(
      xI, sfm_data->GetViews().at(iIndex)->ui_width, sfm_data->GetViews().at(iIndex)->ui_height,
      xJ, sfm_data->GetViews().at(jIndex)->ui_width, sfm_data->GetViews().at(jIndex)->ui_height,
      false); // configure as point to point error model.

    // Robustly estimate the Homography matrix with A Contrario ransac
    const double upper_bound_precision = Square(m_dPrecision);
    std::vector<size_t> vec_inliers;
    const std::pair<double,double> ACRansacOut =
      ACRANSAC(kernel, vec_inliers, m_stIteration, &m_H, upper_bound_precision);

    if (vec_inliers.size() > KernelType::MINIMUM_SAMPLES *2.5)  {
      m_dPrecision_robust = ACRansacOut.first;
      // update geometric_inliers
      geometric_inliers.reserve(vec_inliers.size());
      for ( const size_t & index : vec_inliers)  {
        geometric_inliers.push_back( vec_PutativeMatches[index] );
      }
      return true;
    }
    else  {
      vec_inliers.clear();
      return false;
    }
    return true;
  }

  /// Export point feature based vector to a matrix [(x,y)'T, (x,y)'T]
  /// Use the camera intrinsics in order to get undistorted pixel coordinates
  template<typename MatT >
  static void PointsToMat(
    const cameras::IntrinsicBase * cam,
    const features::PointFeatures & vec_feats,
    MatT & m)
  {
    m.resize(2, vec_feats.size());
    typedef typename MatT::Scalar Scalar; // Output matrix type

    size_t i = 0;
    for( features::PointFeatures::const_iterator iter = vec_feats.begin();
      iter != vec_feats.end(); ++iter, ++i)
    {
      if (cam)
        m.col(i) = cam->get_ud_pixel(Vec2(iter->x(), iter->y()));
      else
        m.col(i) = iter->coords().cast<Scalar>();
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
    if (m_dPrecision_robust != std::numeric_limits<double>::infinity())
    {
      // Get back corresponding view index
      const IndexT iIndex = pairIndex.first;
      const IndexT jIndex = pairIndex.second;

      const sfm::View * view_I = sfm_data->views.at(iIndex).get();
      const sfm::View * view_J = sfm_data->views.at(jIndex).get();

      // Retrieve corresponding pair camera intrinsic if any
      const cameras::IntrinsicBase * cam_I =
        sfm_data->GetIntrinsics().count(view_I->id_intrinsic) ?
          sfm_data->GetIntrinsics().at(view_I->id_intrinsic).get() : nullptr;
      const cameras::IntrinsicBase * cam_J =
        sfm_data->GetIntrinsics().count(view_J->id_intrinsic) ?
          sfm_data->GetIntrinsics().at(view_J->id_intrinsic).get() : nullptr;

      if (dDistanceRatio < 0)
      {
        // Filtering based only on region positions
        const features::PointFeatures pointsFeaturesI = regions_provider->regions_per_view.at(iIndex)->GetRegionsPositions();
        const features::PointFeatures pointsFeaturesJ = regions_provider->regions_per_view.at(jIndex)->GetRegionsPositions();
        Mat xI, xJ;
        PointsToMat(cam_I, pointsFeaturesI, xI);
        PointsToMat(cam_J, pointsFeaturesJ, xJ);

        geometry_aware::GuidedMatching
          <Mat3, openMVG::homography::kernel::AsymmetricError>(
          m_H, xI, xJ, Square(m_dPrecision_robust), matches);

        // Remove duplicates
        matching::IndMatch::getDeduplicated(matches);

        // Remove matches that have the same (X,Y) coordinates
        matching::IndMatchDecorator<float> matchDeduplicator(matches, pointsFeaturesI, pointsFeaturesJ);
        matchDeduplicator.getDeduplicated(matches);
      }
      else
      {
        // Filtering based on region positions and regions descriptors
        geometry_aware::GuidedMatching
          <Mat3, openMVG::homography::kernel::AsymmetricError>(
          m_H,
          cam_I, *regions_provider->regions_per_view.at(iIndex),
          cam_J, *regions_provider->regions_per_view.at(jIndex),
          Square(m_dPrecision_robust), Square(dDistanceRatio),
          matches);
      }
    }
    return matches.size() != 0;
  }

  double m_dPrecision;  //upper_bound precision used for robust estimation
  size_t m_stIteration; //maximal number of iteration for robust estimation
  //
  //-- Stored data
  Mat3 m_H;
  double m_dPrecision_robust;
};

} //namespace matching_image_collection 
} // namespace openMVG


