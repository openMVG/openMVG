
// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#pragma once

#include "openMVG/features/features.hpp"
#include "openMVG/matching/indMatchDecoratorXY.hpp"
#include "openMVG/matching/matching_filters.hpp"
#include "openMVG/matching_image_collection/Matcher.hpp"

#include "nonFree/matching_cascade_hashing/CasHash.hpp"

#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "third_party/progress/progress.hpp"

// [1] Fast and Accurate Image Matching with Cascade Hashing for 3D Reconstruction.
// Authors: Jian Cheng, Cong Leng, Jiaxiang Wu, Hainan Cui, Hanqing Lu.
// Conference: CVPR 2014.

namespace openMVG {

using namespace openMVG::matching;

/// Implementation of an Image Collection Matcher
/// - Cascade Hashing matching [1]
///
template <typename KeypointSetT>
class Matcher_CascadeHashing_AllInMemory
{
  // Alias to internal stored Feature and Descriptor type
  typedef typename KeypointSetT::FeatureT FeatureT;
  typedef typename KeypointSetT::DescriptorT DescriptorT;
  typedef std::vector<DescriptorT > DescsT; // A collection of descriptors
  // Alias to Descriptor value type
  typedef typename DescriptorT::bin_type DescBin_typeT;

  public:
  Matcher_CascadeHashing_AllInMemory(float distRatio):
    fDistRatio(distRatio)
  {
  }

  /// Load all features and descriptors in memory
  bool loadData(
    const std::vector<std::string> & vec_fileNames, // input filenames
    const std::string & sMatchDir) // where the data are saved
  {
    bool bOk = true;
    for (size_t j = 0; j < vec_fileNames.size(); ++j)  {
      // Load descriptor of Jnth image
      const std::string sFeatJ = stlplus::create_filespec(sMatchDir,
        stlplus::basename_part(vec_fileNames[j]), "feat");
      const std::string sDescJ = stlplus::create_filespec(sMatchDir,
        stlplus::basename_part(vec_fileNames[j]), "desc");

      bOk &= loadFeatsFromFile(sFeatJ, map_Feat[j]);
      bOk &= loadDescsFromBinFile(sDescJ, map_Desc[j]);
    }
    if (bOk)
    {
      nonFree::CASHASH::ImportFeatures(map_Desc, vec_hashing);
    }
    return bOk;
  }

  void Match(
    const std::vector<std::string> & vec_fileNames, // input filenames,
    const PairsT & pairs,
    PairWiseMatches & map_PutativesMatches) // the pairwise photometric corresponding points
  {
    C_Progress_display my_progress_bar( pairs.size() );

    // Sort pairs according the first index to minimize memory exchange
    std::map<size_t, std::vector<size_t> > map_Pairs;
    for (PairsT::const_iterator iter = pairs.begin(); iter != pairs.end(); ++iter)
    {
      map_Pairs[iter->first].push_back(iter->second);
    }

    for (std::map<size_t, std::vector<size_t> >::const_iterator iter = map_Pairs.begin();
      iter != map_Pairs.end(); ++iter)
    {
      const size_t I = iter->first;
      const std::vector<FeatureT> & featureSetI = map_Feat[I];

      const std::vector<size_t> & indexToCompare = iter->second;
#ifdef USE_OPENMP
  #pragma omp parallel for schedule(dynamic)
#endif
      for (int j = 0; j < (int)indexToCompare.size(); ++j)
      {
        const size_t J = indexToCompare[j];
        const std::vector<FeatureT> & featureSetJ = map_Feat[J];

        std::vector<IndMatch> vec_FilteredMatches;
        cascadeHashing.MatchSpFast(
          vec_FilteredMatches,
          vec_hashing[I], map_Desc[I],
          vec_hashing[J], map_Desc[J],
          fDistRatio);

        // Remove duplicates
        IndMatch::getDeduplicated(vec_FilteredMatches);

        // Remove matches that have the same X,Y coordinates
        IndMatchDecorator<float> matchDeduplicator(vec_FilteredMatches, featureSetI, featureSetJ);
        matchDeduplicator.getDeduplicated(vec_FilteredMatches);

#ifdef USE_OPENMP
  #pragma omp critical
#endif
        {
          if (!vec_FilteredMatches.empty())
            map_PutativesMatches.insert( make_pair( make_pair(I,J), vec_FilteredMatches ));
          ++my_progress_bar;
        }
      }
    }
  }

  private:
  // Features per image
  std::map<size_t, std::vector<FeatureT> > map_Feat;
  // Descriptors per image as contiguous memory
  std::map<size_t, DescsT > map_Desc;
  // Distance ratio used to discard spurious correspondence
  float fDistRatio;

  // CascadeHashing object
  nonFree::CASHASH::CasHashMatcher cascadeHashing;
  std::vector<nonFree::CASHASH::ImageFeatures> vec_hashing;
};

}; // namespace openMVG
