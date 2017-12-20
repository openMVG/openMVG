// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/matching_image_collection/Geometric_Filter_utils.hpp"

#include "openMVG/cameras/Camera_Intrinsics.hpp"
#include "openMVG/features/feature.hpp"
#include "openMVG/matching/indMatch.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_view.hpp"
#include "openMVG/sfm/pipelines/sfm_features_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_regions_provider.hpp"

namespace openMVG {
namespace matching_image_collection {

void MatchesPointsToMat
(
  const matching::IndMatches & putativeMatches,
  const cameras::IntrinsicBase * cam_I,
  const features::PointFeatures & feature_I,
  const cameras::IntrinsicBase * cam_J,
  const features::PointFeatures & feature_J,
  Mat2X & x_I,
  Mat2X & x_J
)
{
  const size_t n = putativeMatches.size();
  x_I.resize(2, n);
  x_J.resize(2, n);
  using Scalar = typename Mat::Scalar; // Output matrix type

  for (size_t i=0; i < putativeMatches.size(); ++i)  {
    const features::PointFeature & pt_I = feature_I[putativeMatches[i].i_];
    const features::PointFeature & pt_J = feature_J[putativeMatches[i].j_];
    if (cam_I)
      x_I.col(i) = cam_I->get_ud_pixel(pt_I.coords().cast<Scalar>());
    else
      x_I.col(i) = pt_I.coords().cast<Scalar>();

    if (cam_J)
      x_J.col(i) = cam_J->get_ud_pixel(pt_J.coords().cast<Scalar>());
    else
      x_J.col(i) = pt_J.coords().cast<Scalar>();
  }
}

void MatchesPairToMat
(
  const Pair pairIndex,
  const matching::IndMatches & putativeMatches,
  const sfm::SfM_Data * sfm_data,
  const std::shared_ptr<sfm::Regions_Provider> & regions_provider,
  Mat2X & x_I,
  Mat2X & x_J
)
{
  const sfm::View
    * view_I = sfm_data->views.at(pairIndex.first).get(),
    * view_J = sfm_data->views.at(pairIndex.second).get();

  // Retrieve corresponding pair camera intrinsic if any
  const cameras::IntrinsicBase
    * cam_I =
      sfm_data->GetIntrinsics().count(view_I->id_intrinsic) ?
        sfm_data->GetIntrinsics().at(view_I->id_intrinsic).get() : nullptr,
    * cam_J =
      sfm_data->GetIntrinsics().count(view_J->id_intrinsic) ?
        sfm_data->GetIntrinsics().at(view_J->id_intrinsic).get() : nullptr;

  // Load features of Inth and Jnth images
  const std::shared_ptr<features::Regions>
    regionsI = regions_provider->get(pairIndex.first),
    regionsJ = regions_provider->get(pairIndex.second);
  const features::PointFeatures
    feature_I = regionsI->GetRegionsPositions(),
    feature_J = regionsJ->GetRegionsPositions();

  MatchesPointsToMat(
    putativeMatches,
    cam_I, feature_I,
    cam_J, feature_J,
    x_I, x_J);
}

void MatchesPairToMat
(
  const Pair pairIndex,
  const matching::IndMatches & putativeMatches,
  const sfm::SfM_Data * sfm_data,
  const std::shared_ptr<sfm::Features_Provider> & features_provider,
  Mat2X & x_I,
  Mat2X & x_J
)
{
  const sfm::View
    * view_I = sfm_data->views.at(pairIndex.first).get(),
    * view_J = sfm_data->views.at(pairIndex.second).get();

  // Retrieve corresponding pair camera intrinsic if any
  const cameras::IntrinsicBase
    * cam_I =
      sfm_data->GetIntrinsics().count(view_I->id_intrinsic) ?
        sfm_data->GetIntrinsics().at(view_I->id_intrinsic).get() : nullptr,
    * cam_J =
      sfm_data->GetIntrinsics().count(view_J->id_intrinsic) ?
        sfm_data->GetIntrinsics().at(view_J->id_intrinsic).get() : nullptr;

  // Load features of Inth and Jnth images
  const features::PointFeatures
    feature_I = features_provider->feats_per_view.at(pairIndex.first),
    feature_J = features_provider->feats_per_view.at(pairIndex.second);

  MatchesPointsToMat(
    putativeMatches,
    cam_I, feature_I,
    cam_J, feature_J,
    x_I, x_J);
}

} //namespace matching_image_collection
} // namespace openMVG
