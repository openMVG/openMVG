// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MATCHING_IMAGE_COLLECTION_GEOMETRIC_FILTER_HPP
#define OPENMVG_MATCHING_IMAGE_COLLECTION_GEOMETRIC_FILTER_HPP

#include <algorithm>
#include <map>
#include <vector>

#include "openMVG/features/feature.hpp"
#include "openMVG/matching/indMatch.hpp"

#include "third_party/progress/progress_display.hpp"

namespace openMVG { namespace sfm { struct Regions_Provider; } }
namespace openMVG { namespace sfm { struct SfM_Data; } }

namespace openMVG {

namespace matching_image_collection {

using namespace openMVG::matching;

/// Allow to keep only geometrically coherent matches
/// -> It discards pairs that do not lead to a valid robust model estimation
struct ImageCollectionGeometricFilter
{
  ImageCollectionGeometricFilter
  (
    const sfm::SfM_Data * sfm_data,
    const std::shared_ptr<sfm::Regions_Provider> & regions_provider
  ):sfm_data_(sfm_data),
    regions_provider_(regions_provider)
  {}

  /// Perform robust model estimation (with optional guided_matching) for all
  /// the pairs and regions correspondences contained in the putative_matches
  /// set.
  template<typename GeometryFunctor>
  void Robust_model_estimation
  (
    const GeometryFunctor & functor,
    const PairWiseMatches & putative_matches,
    const bool b_guided_matching = false,
    const double d_distance_ratio = 0.6,
    C_Progress *progress_bar = nullptr
  );

  const PairWiseMatches & Get_geometric_matches() const
  {
    return _map_GeometricMatches;
  }

  // Data
  const sfm::SfM_Data * sfm_data_;
  const std::shared_ptr<sfm::Regions_Provider> & regions_provider_;
  PairWiseMatches _map_GeometricMatches;
};

template<typename GeometryFunctor>
void ImageCollectionGeometricFilter::Robust_model_estimation
(
  const GeometryFunctor & functor,
  const PairWiseMatches & putative_matches,
  const bool b_guided_matching,
  const double d_distance_ratio,
  C_Progress * my_progress_bar
)
{
  if (!my_progress_bar)
    my_progress_bar = &C_Progress::dummy();
  my_progress_bar->restart( putative_matches.size(), "\n- Geometric filtering -\n" );

#ifdef OPENMVG_USE_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (int i = 0; i < (int)putative_matches.size(); ++i)
  {
    if (my_progress_bar->hasBeenCanceled())
      continue;
    auto iter = putative_matches.begin();
    advance(iter,i);

    Pair current_pair = iter->first;
    const std::vector<IndMatch> & vec_PutativeMatches = iter->second;

    //-- Apply the geometric filter (robust model estimation)
    {
      IndMatches putative_inliers;
      GeometryFunctor geometricFilter = functor; // use a copy since we are in a multi-thread context
      if (geometricFilter.Robust_estimation(
        sfm_data_,
        regions_provider_,
        iter->first,
        vec_PutativeMatches,
        putative_inliers))
      {
        if (b_guided_matching)
        {
          IndMatches guided_geometric_inliers;
          geometricFilter.Geometry_guided_matching(
            sfm_data_,
            regions_provider_,
            iter->first,
            d_distance_ratio,
            guided_geometric_inliers);
          //std::cout
          // << "#before/#after: " << putative_inliers.size()
          // << "/" << guided_geometric_inliers.size() << std::endl;
          std::swap(putative_inliers, guided_geometric_inliers);
        }

#ifdef OPENMVG_USE_OPENMP
#pragma omp critical
#endif
        {
          _map_GeometricMatches.insert( {current_pair, std::move(putative_inliers)});
        }
      }
    }
    ++(*my_progress_bar);
  }
}

} // namespace matching_image_collection
} // namespace openMVG

#endif // OPENMVG_MATCHING_IMAGE_COLLECTION_GEOMETRIC_FILTER_HPP
