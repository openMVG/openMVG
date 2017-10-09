// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/matching_image_collection/Matcher_Regions.hpp"
#include "openMVG/matching_image_collection/Matcher.hpp"
#include "openMVG/matching/regions_matcher.hpp"
#include "openMVG/sfm/pipelines/sfm_regions_provider.hpp"

#include "third_party/progress/progress.hpp"

namespace openMVG {
namespace matching_image_collection {

using namespace openMVG::matching;
using namespace openMVG::features;

void Matcher_Regions::Match(
  const sfm::SfM_Data & sfm_data,
  const std::shared_ptr<sfm::Regions_Provider> & regions_provider,
  const Pair_Set & pairs,
  PairWiseMatchesContainer & map_PutativesMatches,
  C_Progress * my_progress_bar)const
{
  if (!create_matcher_ || !call_matcher_)
    return;

  if (!my_progress_bar)
    my_progress_bar = &C_Progress::dummy();
#ifdef OPENMVG_USE_OPENMP
  std::cout << "Using the OPENMP thread interface" << std::endl;
  //const bool b_multithreaded_pair_search = (eMatcherType_ == CASCADE_HASHING_L2);
  // -> set to true for CASCADE_HASHING_L2, since OpenMP instructions are not used in this matcher
#endif

  my_progress_bar->restart(pairs.size(), "\n- Matching -\n");

  // Sort pairs according the first index to minimize the MatcherT build operations
  using Map_vectorT = std::map<IndexT, std::vector<IndexT>>;
  Map_vectorT map_Pairs;
  for (const auto & pair_idx : pairs)
  {
    map_Pairs[pair_idx.first].push_back(pair_idx.second);
  }

  // Perform matching between all the pairs
  for (const auto & pairs : map_Pairs)
  {
    if (my_progress_bar->hasBeenCanceled())
      continue;
    const IndexT I = pairs.first;
    const auto & indexToCompare = pairs.second;

    const auto regionsI = regions_provider->get(I);
    if (regionsI->RegionCount() == 0)
    {
      (*my_progress_bar) += indexToCompare.size();
      continue;
    }

    // Initialize the matching interface
    const auto matcher = create_matcher_(*regionsI.get());
    if (!matcher)
      continue;

#ifdef OPENMVG_USE_OPENMP
    //#pragma omp parallel for schedule(dynamic) if (b_multithreaded_pair_search)
#endif
    for (int j = 0; j < (int)indexToCompare.size(); ++j)
    {
      const IndexT J = indexToCompare[j];

      const auto regionsJ = regions_provider->get(J);
      if (regionsJ->RegionCount() == 0
          || regionsI->Type_id() != regionsJ->Type_id())
      {
        ++(*my_progress_bar);
        continue;
      }

      IndMatches matches;
      call_matcher_(matcher.get(), *regionsJ.get(), matches);

      if (call_post_process_match_functor)
      {
        matches = call_post_process_match_functor(
          *regionsI.get(),
          *regionsJ.get(),
          matches,
          {sfm_data.GetViews().at(I)->ui_width, sfm_data.GetViews().at(I)->ui_height},
          {sfm_data.GetViews().at(J)->ui_width, sfm_data.GetViews().at(J)->ui_height});
      }

#ifdef OPENMVG_USE_OPENMP
  #pragma omp critical
#endif
      {
        if (!matches.empty())
        {
          map_PutativesMatches.insert( { {I,J}, std::move(matches) } );
        }
      }
      ++(*my_progress_bar);
    }
  }
}

} // namespace matching_image_collection
} // namespace openMVG
