
// Copyright (c) 2012, 2013, 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#pragma once

#include "openMVG/features/features.hpp"
#include "openMVG/matching/indMatchDecoratorXY.hpp"
#include "openMVG/matching/matching_filters.hpp"
#include "openMVG/matching_image_collection/Matcher.hpp"

#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "third_party/progress/progress.hpp"

namespace openMVG {

using namespace openMVG::matching;

/// Implementation of an Image Collection Matcher
/// Compute putative matches between a collection of pictures
/// Spurious correspondences are discarded by using the
///  a threshold over the distance ratio of the 2 neighbours points.
///
template <typename KeypointSetT, typename MatcherT>
class Matcher_AllInMemory : public Matcher
{
  // Alias to internal stored Feature and Descriptor type
  typedef typename KeypointSetT::FeatureT FeatureT;
  typedef typename KeypointSetT::DescriptorT DescriptorT;
  typedef std::vector<DescriptorT > DescsT; // A collection of descriptors
  // Alias to Descriptor value type
  typedef typename DescriptorT::bin_type DescBin_typeT;

  public:
  Matcher_AllInMemory(float distRatio):
    Matcher(),
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
    return bOk;
  }

  void Match(
    const std::vector<std::string> & vec_fileNames, // input filenames,
    const PairsT & pairs,
    PairWiseMatches & map_PutativesMatches)const // the pairwise photometric corresponding points
  {
#ifdef OPENMVG_USE_OPENMP
    std::cout << "Using the OPENMP thread interface" << std::endl;
#endif
    C_Progress_display my_progress_bar( pairs.size() );

    // Sort pairs according the first index to minimize the MatcherT build operations
    std::map<size_t, std::vector<size_t> > map_Pairs;
    for (PairsT::const_iterator iter = pairs.begin(); iter != pairs.end(); ++iter)
    {
      map_Pairs[iter->first].push_back(iter->second);
    }

    for (std::map<size_t, std::vector<size_t> >::const_iterator iter = map_Pairs.begin();
      iter != map_Pairs.end(); ++iter)
    {
      const size_t I = iter->first;
      // Load features and descriptors of Inth image
      typename std::map<size_t, std::vector<FeatureT> >::const_iterator iter_FeaturesI = map_Feat.find(I);
      typename std::map<size_t, DescsT >::const_iterator iter_DescriptorI = map_Desc.find(I);

      const std::vector<FeatureT> & featureSetI = iter_FeaturesI->second;
      const size_t featureSetI_Size = iter_FeaturesI->second.size();
      const DescBin_typeT * tab0 =
        reinterpret_cast<const DescBin_typeT *>(&iter_DescriptorI->second[0]);

      MatcherT matcher10;
      ( matcher10.Build(tab0, featureSetI_Size, DescriptorT::static_size) );

      const std::vector<size_t> & indexToCompare = iter->second;
#ifdef OPENMVG_USE_OPENMP
  #pragma omp parallel for schedule(dynamic)
#endif
      for (int j = 0; j < (int)indexToCompare.size(); ++j)
      {
        const size_t J = indexToCompare[j];
        // Load descriptors of Jnth image
        typename std::map<size_t, std::vector<FeatureT> >::const_iterator iter_FeaturesJ = map_Feat.find(J);
        typename std::map<size_t, DescsT >::const_iterator iter_DescriptorJ = map_Desc.find(J);

        const std::vector<FeatureT> & featureSetJ = iter_FeaturesJ->second;
        const DescBin_typeT * tab1 =
          reinterpret_cast<const DescBin_typeT *>(&iter_DescriptorJ->second[0]);

        const size_t NNN__ = 2;
        std::vector<int> vec_nIndice10;
        std::vector<typename MatcherT::DistanceType> vec_fDistance10;

        //Find left->right
        matcher10.SearchNeighbours(tab1, featureSetJ.size(), &vec_nIndice10, &vec_fDistance10, NNN__);

        std::vector<IndMatch> vec_FilteredMatches;
        std::vector<int> vec_NNRatioIndexes;
        NNdistanceRatio( vec_fDistance10.begin(), // distance start
          vec_fDistance10.end(),  // distance end
          NNN__, // Number of neighbor in iterator sequence (minimum required 2)
          vec_NNRatioIndexes, // output (index that respect Lowe Ratio)
          Square(fDistRatio)); // squared dist ratio due to usage of a squared metric

        for (size_t k=0; k < vec_NNRatioIndexes.size()-1&& vec_NNRatioIndexes.size()>0; ++k)
        {
          const size_t index = vec_NNRatioIndexes[k];
          vec_FilteredMatches.push_back(
            IndMatch(vec_nIndice10[index*NNN__], index) );
        }

        // Remove duplicates
        IndMatch::getDeduplicated(vec_FilteredMatches);

        // Remove matches that have the same X,Y coordinates
        IndMatchDecorator<float> matchDeduplicator(vec_FilteredMatches, featureSetI, featureSetJ);
        matchDeduplicator.getDeduplicated(vec_FilteredMatches);

#ifdef OPENMVG_USE_OPENMP
  #pragma omp critical
#endif
        {
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
};

}; // namespace openMVG
