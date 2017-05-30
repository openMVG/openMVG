// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Matthew Daiter.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/matching_image_collection/GPU_Matcher_Regions.hpp"

#ifdef OPENMVG_USE_CUDA
#include "openMVG/matching_image_collection/gpu/LatchBitMatcher.hpp"
#endif

#include "openMVG/matching_image_collection/Matcher.hpp"
#include "openMVG/matching_image_collection/Matcher_Regions.hpp"
#include "openMVG/matching/regions_matcher.hpp"
#include "openMVG/sfm/pipelines/sfm_regions_provider.hpp"

#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "third_party/progress/progress.hpp"

namespace openMVG {
namespace matching_image_collection {

using namespace openMVG::matching;
using namespace openMVG::features;

GPU_Matcher_Regions::GPU_Matcher_Regions(
  float distRatio, EMatcherType eMatcherType)
  :Matcher(), f_dist_ratio_(distRatio), eMatcherType_(eMatcherType)
{
}

void GPU_Matcher_Regions::Match(
  const sfm::SfM_Data & sfm_data,
  const std::shared_ptr<sfm::Regions_Provider> & regions_provider,
  const Pair_Set & pairs,
  PairWiseMatchesContainer & map_PutativesMatches,
  C_Progress * my_progress_bar)const // the pairwise photometric corresponding points
{
#ifdef OPENMVG_USE_OPENMP
  std::cout << "Using the OPENMP thread interface" << std::endl;
#endif

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

    std::shared_ptr<features::Regions> regionsI = regions_provider->get(I);
    if (regionsI->RegionCount() == 0)
    {
      my_progress_bar += indexToCompare.size();
      continue;
    }

    switch (eMatcherType_) {
      case BRUTE_FORCE_HAMMING:
      {
        std::vector<IndMatches> matchedPoints(indexToCompare.size());
#ifdef OPENMVG_USE_OPENMP
        #pragma omp parallel for schedule(dynamic)
#endif
        for (size_t j = 0; j < (int)indexToCompare.size(); ++j)
        {
          const size_t J = indexToCompare[j];

          std::shared_ptr<features::Regions> regionsJ = regions_provider->get(J);
          if (regionsJ->RegionCount() == 0
              || regionsI->Type_id() != regionsJ->Type_id())
          {
#ifdef OPENMVG_USE_OPENMP
      #pragma omp critical
#endif
            ++my_progress_bar;
            continue;
          }

          switch (regionsI->DescriptorLength()) {
            case 64: {
#ifdef OPENMVG_USE_CUDA
              gpu::LatchBitMatcher matcher(std::max(regionsI->RegionCount(), regionsJ->RegionCount()));
              matcher.match(
                const_cast<unsigned int*>(static_cast<const unsigned int*>(regionsI->DescriptorRawData())),
                const_cast<unsigned int*>(static_cast<const unsigned int*>(regionsJ->DescriptorRawData())), 
                regionsI->RegionCount(), 
                regionsJ->RegionCount());
              matchedPoints[j] = matcher.retrieveMatches(f_dist_ratio_);
#endif
              break;
            }
            default: {
              std::cout << "NOT PERFORMING GPU MATCHING: "
                        << "cannot handle this size of descriptor"
                        << std::endl;
              break;
            }
          }
        }
#ifdef OPENMVG_USE_OPENMP
        #pragma omp parallel for schedule(dynamic)
#endif
        for (int j = 0; j < static_cast<int>(indexToCompare.size()); j++)
        {
          const size_t J = indexToCompare[j];

          IndMatches vec_putatives_matches;
          for (size_t k = 0; k < matchedPoints[j].size(); k++) {
            vec_putatives_matches.emplace_back(matchedPoints[j][k].i_, matchedPoints[j][k].j_);
          }
          ++my_progress_bar;
#ifdef OPENMVG_USE_OPENMP
          #pragma omp critical
#endif
          {
            if (!vec_putatives_matches.empty())
            {
              map_PutativesMatches.insert( { {I,J}, std::move(vec_putatives_matches) });
            }
          }
        }
      }
      default:
      {
        break;
      }
    }
  }
}

} // namespace openMVG
} // namespace matching_image_collection
