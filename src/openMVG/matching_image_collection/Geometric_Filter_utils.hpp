
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

namespace openMVG {
namespace matching_image_collection {

template<typename MatT >
void MatchesPointsToMat
(
  const matching::IndMatches & putativeMatches,
  const cameras::IntrinsicBase * cam_I,
  const features::PointFeatures & feature_I,
  const cameras::IntrinsicBase * cam_J,
  const features::PointFeatures & feature_J,
  MatT & x_I, MatT & x_J
)
{
  const size_t n = putativeMatches.size();
  x_I.resize(2, n);
  x_J.resize(2, n);
  typedef typename MatT::Scalar Scalar; // Output matrix type

  for (size_t i=0; i < putativeMatches.size(); ++i)  {
    const features::PointFeature & pt_I = feature_I[putativeMatches[i]._i];
    const features::PointFeature & pt_J = feature_J[putativeMatches[i]._j];
    if (cam_I)
      x_I.col(i) = cam_I->get_ud_pixel(pt_I.coords().cast<double>());
    else
      x_I.col(i) = pt_I.coords().cast<double>();

    if (cam_J)
      x_J.col(i) = cam_J->get_ud_pixel(pt_J.coords().cast<double>());
    else
      x_J.col(i) = pt_J.coords().cast<double>();
  }
}

template<typename MatT >
void MatchesPairToMat
(
  const Pair pairIndex,
  const matching::IndMatches & putativeMatches,
  const sfm::SfM_Data * sfm_data,
  const std::shared_ptr<sfm::Regions_Provider> & regions_provider,
  MatT & x_I, MatT & x_J
)
{
  const sfm::View * view_I = sfm_data->views.at(pairIndex.first).get();
  const sfm::View * view_J = sfm_data->views.at(pairIndex.second).get();

  // Retrieve corresponding pair camera intrinsic if any
  const cameras::IntrinsicBase * cam_I =
    sfm_data->GetIntrinsics().count(view_I->id_intrinsic) ?
      sfm_data->GetIntrinsics().at(view_I->id_intrinsic).get() : NULL;
  const cameras::IntrinsicBase * cam_J =
    sfm_data->GetIntrinsics().count(view_J->id_intrinsic) ?
      sfm_data->GetIntrinsics().at(view_J->id_intrinsic).get() : NULL;

  // Load features of Inth and Jnth images
  const features::PointFeatures feature_I = regions_provider->regions_per_view.at(pairIndex.first)->GetRegionsPositions();
  const features::PointFeatures feature_J = regions_provider->regions_per_view.at(pairIndex.second)->GetRegionsPositions();

  MatchesPointsToMat(
    putativeMatches,
    cam_I, feature_I,
    cam_J, feature_J,
    x_I, x_J);
}

} // namespace openMVG
} //namespace matching_image_collection
