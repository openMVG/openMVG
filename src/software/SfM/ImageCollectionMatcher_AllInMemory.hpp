
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#pragma once

#include "software/SfM/ImageCollectionMatcher.hpp"
#include "openMVG/features/features.hpp"
#include "openMVG/matching/indMatchDecoratorXY.hpp"
#include "openMVG/matching/matching_filters.hpp"


using namespace openMVG;

#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "third_party/progress/progress.hpp"

/// Implementation of an Image Collection Matcher
/// Compute putative matches between a collection of pictures
/// Suppose a symmetric matching results : Will compare the
///  upper matrix index for image matching.
/// For descriptor matching between two image indexes :
///  The distance ratio of the 2 neighbours points is used
///   to discard spurious correspondences.
template <typename KeypointSetT, typename MatcherT>
class ImageCollectionMatcher_AllInMemory : public ImageCollectionMatcher
{
  // Alias to internal stored Feature and Descriptor type
  typedef typename KeypointSetT::FeatureT FeatureT;
  typedef typename KeypointSetT::DescriptorT DescriptorT;
  typedef std::vector<DescriptorT > DescsT; // A collection of descriptors
  // Alias to Descriptor value type
  typedef typename DescriptorT::bin_type DescBin_typeT;

  public:
  ImageCollectionMatcher_AllInMemory(float distRatio) :ImageCollectionMatcher(), fDistRatio(distRatio)
  {
  }

  /// Load all features and descriptors in memory
  bool loadData(
    const std::vector<std::string> & vec_fileNames, // input filenames
    const std::string & sMatchDir) // where the data are saved
  {
    // Alias for descriptor internal data type and length

    const int static_size = DescriptorT::static_size;

    bool bOk = true;
    for (size_t j = 0; j < vec_fileNames.size(); ++j)  {
      // Load descriptor of Jnth image
      const std::string sFeatJ = stlplus::create_filespec(sMatchDir,
        stlplus::basename_part(vec_fileNames[j]), "feat");

      const std::string sDescJ = stlplus::create_filespec(sMatchDir,
        stlplus::basename_part(vec_fileNames[j]), "desc");
      bOk &= loadFeatsFromFile(sFeatJ, map_Feat[j]);

      KeypointSetT kpSetI;
      bOk &= loadDescsFromBinFile(sDescJ, map_Desc[j]);
    }
    return bOk;
  }

  void Match(
    const std::vector<std::string> & vec_fileNames, // input filenames,
    IndexedMatchPerPair & map_PutativesMatches)const // the pairwise photometric corresponding points
  {
#ifdef USE_OPENMP
    std::cout << "Using the OPENMP thread interface" << std::endl;
#endif
    C_Progress_display my_progress_bar( vec_fileNames.size()*(vec_fileNames.size()-1) / 2.0 );

    for (size_t i = 0; i < vec_fileNames.size(); ++i)
    {
      // Load features and descriptors of Inth image
      typename std::map<size_t, std::vector<FeatureT> >::const_iterator iter_FeaturesI = map_Feat.begin();
      typename std::map<size_t, DescsT >::const_iterator iter_DescriptorI = map_Desc.begin();
      std::advance(iter_FeaturesI, i);
      std::advance(iter_DescriptorI, i);

      const std::vector<FeatureT> & featureSetI = iter_FeaturesI->second;
      const size_t featureSetI_Size = iter_FeaturesI->second.size();
      const DescBin_typeT * tab0 =
        reinterpret_cast<const DescBin_typeT *>(&iter_DescriptorI->second[0]);

      MatcherT matcher10;
      ( matcher10.Build(tab0, featureSetI_Size, DescriptorT::static_size) );

#ifdef USE_OPENMP
  #pragma omp parallel for schedule(dynamic, 1)
#endif
      for (int j = i+1; j < (int)vec_fileNames.size(); ++j)
      {
        // Load descriptor of Jnth image
        typename std::map<size_t, std::vector<FeatureT> >::const_iterator iter_FeaturesJ = map_Feat.begin();
        typename std::map<size_t, DescsT >::const_iterator iter_DescriptorJ = map_Desc.begin();
        std::advance(iter_FeaturesJ, j);
        std::advance(iter_DescriptorJ, j);

        const std::vector<FeatureT> & featureSetJ = iter_FeaturesJ->second;
        const size_t featureSetJ_Size = iter_FeaturesJ->second.size();
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
          vec_FilteredMatches.push_back(
            IndMatch(vec_nIndice10[vec_NNRatioIndexes[k]*NNN__],
                     vec_NNRatioIndexes[k]) );
        }

        // Remove duplicates
        IndMatch::getDeduplicated(vec_FilteredMatches);

        // Remove matches that have the same X,Y coordinates
        IndMatchDecorator<float> matchDeduplicator(vec_FilteredMatches, featureSetI, featureSetJ);
        matchDeduplicator.getDeduplicated(vec_FilteredMatches);

#ifdef USE_OPENMP
  #pragma omp critical
#endif
        {
          map_PutativesMatches.insert( make_pair( make_pair(i,j), vec_FilteredMatches ));
        }

        ++my_progress_bar;
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

