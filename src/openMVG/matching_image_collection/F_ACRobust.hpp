
// Copyright (c) 2012, 2013, 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "openMVG/multiview/solver_fundamental_kernel.hpp"
#include "openMVG/multiview/essential.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansac.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansacKernelAdaptator.hpp"

#include "openMVG/matching/indMatch.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/pipelines/sfm_regions_provider.hpp"
#include "openMVG/matching_image_collection/Geometric_Filter_utils.hpp"


namespace openMVG {
namespace matching_image_collection {

//-- A contrario fundamental matrix estimation template functor used for filter pair of putative correspondences 
struct GeometricFilter_FMatrix_AC
{
  GeometricFilter_FMatrix_AC(
    double dPrecision = std::numeric_limits<double>::infinity(),
    size_t iteration = 4096)
    : m_dPrecision(dPrecision), m_stIteration(iteration), m_F(Mat3::Identity()),
      m_dPrecision_robust(std::numeric_limits<double>::infinity()){};

  /// Robust fitting of the FUNDAMENTAL matrix
  bool Robust_estimation(
    const sfm::SfM_Data * sfm_data,
    const std::shared_ptr<sfm::Regions_Provider> & regions_provider,
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

    // Define the AContrario adapted Fundamental matrix solver
    typedef ACKernelAdaptor<
      openMVG::fundamental::kernel::SevenPointSolver,
      openMVG::fundamental::kernel::SimpleError,
      //openMVG::fundamental::kernel::SymmetricEpipolarDistanceError,
      UnnormalizerT,
      Mat3>
      KernelType;

    const KernelType kernel(
      xI, sfm_data->GetViews().at(iIndex)->ui_width, sfm_data->GetViews().at(iIndex)->ui_height,
      xJ, sfm_data->GetViews().at(jIndex)->ui_width, sfm_data->GetViews().at(jIndex)->ui_height, true);

    // Robustly estimate the Fundamental matrix with A Contrario ransac
    const double upper_bound_precision = Square(m_dPrecision);
    std::vector<size_t> vec_inliers;
    const std::pair<double,double> ACRansacOut =
      ACRANSAC(kernel, vec_inliers, m_stIteration, &m_F, upper_bound_precision);

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
  }

  double m_dPrecision;  //upper_bound precision used for robust estimation
  size_t m_stIteration; //maximal number of iteration for robust estimation
  //
  //-- Stored data
  Mat3 m_F;
  double m_dPrecision_robust;
};

} // namespace openMVG
} //namespace matching_image_collection

