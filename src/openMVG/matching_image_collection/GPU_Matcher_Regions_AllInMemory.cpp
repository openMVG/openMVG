
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/matching_image_collection/GPU_Matcher_Regions_AllInMemory.hpp"
#include "openMVG/matching/regions_matcher.hpp"
#include "openMVG/matching_image_collection/Matcher.hpp"

#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "third_party/progress/progress.hpp"

#include "gpu/LatchBitMatcher.hpp"

namespace openMVG {
namespace matching_image_collection {

using namespace openMVG::matching;
using namespace openMVG::features;

GPU_Matcher_Regions_AllInMemory::GPU_Matcher_Regions_AllInMemory(
  float distRatio)
  :Matcher(), f_dist_ratio_(distRatio)
{
}

void GPU_Matcher_Regions_AllInMemory::Match(
  const sfm::SfM_Data & sfm_data,
  const std::shared_ptr<sfm::Regions_Provider> & regions_provider,
  const Pair_Set & pairs,
  PairWiseMatches & map_PutativesMatches)const // the pairwise photometric corresponding points
{
#ifdef OPENMVG_USE_OPENMP
  std::cout << "Using the OPENMP thread interface" << std::endl;
#endif
  C_Progress_display my_progress_bar( pairs.size() );

  // Sort pairs according the first index to minimize the MatcherT build operations
  typedef std::map<size_t, std::vector<size_t> > Map_vectorT;
  Map_vectorT map_Pairs;
  for (Pair_Set::const_iterator iter = pairs.begin(); iter != pairs.end(); ++iter)
  {
    map_Pairs[iter->first].push_back(iter->second);
  }

  // Perform matching between all the pairs
  for (Map_vectorT::const_iterator iter = map_Pairs.begin();
    iter != map_Pairs.end(); ++iter)
  {
    const size_t I = iter->first;
    const std::vector<size_t> & indexToCompare = iter->second;

    const features::Regions & regionsI = *regions_provider->regions_per_view.at(I).get();
    if (regionsI.RegionCount() == 0)
    {
      my_progress_bar += indexToCompare.size();
      continue;
    }
    
    // LatchClassifier for the GPU
    LatchBitMatcher matcher;

#ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for schedule(dynamic)
#endif
    for (int j = 0; j < (int)indexToCompare.size(); ++j)
    {
      const size_t J = indexToCompare[j];

      const features::Regions &regionsJ = *regions_provider->regions_per_view.at(J).get();
      if (regionsJ.RegionCount() == 0
          || regionsI.Type_id() != regionsJ.Type_id())
      {
#ifdef OPENMVG_USE_OPENMP
  #pragma omp critical
#endif
        ++my_progress_bar;
        continue;
      }

      IndMatches vec_putatives_matches;
      // Perform matching in-class. Don't write needless code
   
      auto matchedPoints = matcher.match(
        const_cast<unsigned int*>(static_cast<const unsigned int*>(regionsI.DescriptorRawData())),
        const_cast<unsigned int*>(static_cast<const unsigned int*>(regionsJ.DescriptorRawData())), 
        regionsI.RegionCount(), 
        regionsJ.RegionCount());
      for (size_t k = 0; k < matchedPoints.size(); k++) {
          std::cout << "Matched indices: " << matchedPoints[k].queryIdx << " and " << matchedPoints[k].trainIdx << std::endl;
          vec_putatives_matches.push_back(IndMatch(matchedPoints[k].queryIdx, matchedPoints[k].trainIdx));
      }

#ifdef OPENMVG_USE_OPENMP
  #pragma omp critical
#endif
      {
        ++my_progress_bar;
        if (!vec_putatives_matches.empty())
        {
          map_PutativesMatches.insert( make_pair( make_pair(I,J), std::move(vec_putatives_matches) ));
        }
      }
    }
  }
}

} // namespace openMVG
} // namespace matching_image_collection
