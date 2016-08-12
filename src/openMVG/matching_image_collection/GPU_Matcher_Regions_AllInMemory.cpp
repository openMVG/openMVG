
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
#include "gpu/BruteForceL2Matcher.hpp"

namespace openMVG {
namespace matching_image_collection {

using namespace openMVG::matching;
using namespace openMVG::features;

GPU_Matcher_Regions_AllInMemory::GPU_Matcher_Regions_AllInMemory(
  float distRatio, EMatcherType eMatcherType)
  :Matcher(), f_dist_ratio_(distRatio), eMatcherType_(eMatcherType)
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

		switch (eMatcherType_) {
			case BRUTE_FORCE_HAMMING:
			{
				LatchBitMatcher matchers[indexToCompare.size()];
#ifdef OPENMVG_USE_OPENMP
				omp_set_num_threads(12);
				#pragma omp parallel for schedule(dynamic)
#endif
				for (int j = 0; j < static_cast<int>(indexToCompare.size()); ++j)
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
					// LatchClassifier for the GPU
					matchers[j].match(
						const_cast<unsigned int*>(static_cast<const unsigned int*>(regionsI.DescriptorRawData())),
						const_cast<unsigned int*>(static_cast<const unsigned int*>(regionsJ.DescriptorRawData())), 
						regionsI.RegionCount(), 
						regionsJ.RegionCount());
				}
#ifdef OPENMVG_USE_OPENMP
				#pragma omp parallel for schedule(dynamic)
#endif
				for (int j = 0; j < static_cast<int>(indexToCompare.size()); j++)
				{
					const size_t J = indexToCompare[j];
		 
					std::vector<LatchBitMatcherMatch> matchedPoints = matchers[j].retrieveMatches();
					IndMatches vec_putatives_matches;
					for (size_t k = 0; k < matchedPoints.size(); k++) {
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
				break;
			}
			case BRUTE_FORCE_L2:
			{
			 	std::vector<LatchBitMatcherMatch> matchedPoints[indexToCompare.size()];
#ifdef OPENMVG_USE_OPENMP
				omp_set_num_threads(1);
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
					// LatchClassifier for the GPU
					if (regionsI.Type_id() == typeid(unsigned char).name()) {
						switch (regionsI.DescriptorLength()) {
							case 128: {
								GPUBruteForceL2Matcher<unsigned char, 128> matcher(f_dist_ratio_);
								matcher.match(
									const_cast<unsigned char*>(static_cast<const unsigned char*>(regionsI.DescriptorRawData())),
									const_cast<unsigned char*>(static_cast<const unsigned char*>(regionsJ.DescriptorRawData())), 
									regionsI.RegionCount(), 
									regionsJ.RegionCount());
								matchedPoints[j] = matcher.retrieveMatches();
								break;
							}
							default:
								break;
						}
					}
					else if (regionsI.Type_id() == typeid(float).name()) {
						switch (regionsI.DescriptorLength()) {
							case 128: {
								GPUBruteForceL2Matcher<float, 128> matcher(f_dist_ratio_);
								matcher.match(
									const_cast<float*>(static_cast<const float*>(regionsI.DescriptorRawData())),
									const_cast<float*>(static_cast<const float*>(regionsJ.DescriptorRawData())), 
									regionsI.RegionCount(), 
									regionsJ.RegionCount());
								matchedPoints[j] = matcher.retrieveMatches();
								break;
							}
							case 256: {
								GPUBruteForceL2Matcher<float, 256> matcher(f_dist_ratio_);
								matcher.match(
									const_cast<float*>(static_cast<const float*>(regionsI.DescriptorRawData())),
									const_cast<float*>(static_cast<const float*>(regionsJ.DescriptorRawData())), 
									regionsI.RegionCount(), 
									regionsJ.RegionCount());
								matchedPoints[j] = matcher.retrieveMatches();
								break;
							}
							case 512: {
								GPUBruteForceL2Matcher<float, 512> matcher(f_dist_ratio_);
								matcher.match(
									const_cast<float*>(static_cast<const float*>(regionsI.DescriptorRawData())),
									const_cast<float*>(static_cast<const float*>(regionsJ.DescriptorRawData())), 
									regionsI.RegionCount(), 
									regionsJ.RegionCount());
								matchedPoints[j] = matcher.retrieveMatches();
								break;
							}
							default:
								break;
						}
					}
				}
#ifdef OPENMVG_USE_OPENMP
				#pragma omp parallel for schedule(dynamic)
#endif
				for (int j = 0; j < (int)indexToCompare.size(); j++)
				{
					const size_t J = indexToCompare[j];
		 
					IndMatches vec_putatives_matches;
					for (size_t k = 0; k < matchedPoints[j].size(); k++) {
						vec_putatives_matches.push_back(IndMatch(matchedPoints[j][k].queryIdx, matchedPoints[j][k].trainIdx));
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
				break;
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
