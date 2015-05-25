
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/matching_image_collection/Matcher_Regions_AllInMemory.hpp"

#include "openMVG/features/regions.hpp"
#include "openMVG/matching_image_collection/Matcher.hpp"
#include "openMVG/matching/matcher_brute_force.hpp"
#include "openMVG/matching/matcher_kdtree_flann.hpp"
#include "openMVG/matching/indMatchDecoratorXY.hpp"
#include "openMVG/matching/matching_filters.hpp"

#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "third_party/progress/progress.hpp"

namespace openMVG {

using namespace openMVG::matching;
using namespace openMVG::features;

Matcher_Regions_AllInMemory::Matcher_Regions_AllInMemory(
  float distRatio, EMatcherType eMatcherType)
  :Matcher(), fDistRatio(distRatio), _eMatcherType(eMatcherType)
{
}

/// Load all features and descriptors in memory
bool Matcher_Regions_AllInMemory::loadData(
  std::unique_ptr<features::Regions>& region_type, // interface to load computed regions
  const std::vector<std::string> & vec_fileNames, // input filenames
  const std::string & sMatchDir) // where the data are saved
{
  bool bOk = true;
  for (size_t j = 0; j < vec_fileNames.size(); ++j)  {
    // Load regions of Jnth image
    const std::string sFeatJ = stlplus::create_filespec(sMatchDir,
      stlplus::basename_part(vec_fileNames[j]), "feat");
    const std::string sDescJ = stlplus::create_filespec(sMatchDir,
      stlplus::basename_part(vec_fileNames[j]), "desc");

    regions_perImage[j] = std::unique_ptr<features::Regions>(region_type->EmptyClone());
    bOk &= regions_perImage[j]->Load(sFeatJ, sDescJ);
  }
  return bOk;
}

/// Template matching
template <typename MatcherT>
void Template_Matcher(
  const std::vector<std::string> & vec_fileNames, // input filenames,
  const Pair_Set & pairs,
  PairWiseMatches & map_PutativesMatches, // the pairwise photometric corresponding points
  // Data & parameters
  const std::map<IndexT, std::unique_ptr<features::Regions> > & regions_perImage,
  // Distance ratio used to discard spurious correspondence
  const float fDistRatio,
  // Matcher Type
  EMatcherType _eMatcherType)
{
#ifdef OPENMVG_USE_OPENMP
  std::cout << "Using the OPENMP thread interface" << std::endl;
#endif
  switch(_eMatcherType)
  {
    case BRUTE_FORCE_L2:
      std::cout << "Using BRUTE_FORCE_L2 matcher" << std::endl;
    break;
    case BRUTE_FORCE_HAMMING:
      std::cout << "Using BRUTE_FORCE_HAMMING matcher" << std::endl;
    break;
    case ANN_L2:
      std::cout << "Using ANN_L2 matcher" << std::endl;
    break;
    default:
      std::cout << "Using unknown matcher type" << std::endl;
  }

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

    const features::Regions *regionsI = regions_perImage.at(I).get();
    const size_t regions_countI = regionsI->RegionCount();
    const std::vector<PointFeature> pointFeaturesI = regionsI->GetRegionsPositions();
    const typename MatcherT::ScalarT * tabI = NULL;
    if(regions_countI > 0)
      tabI = reinterpret_cast<const typename MatcherT::ScalarT *>(regionsI->DescriptorRawData());

    MatcherT matcher10;
    ( matcher10.Build(tabI, regions_countI, regionsI->DescriptorLength()) );

    const std::vector<size_t> & indexToCompare = iter->second;

    for (int j = 0; j < (int)indexToCompare.size(); ++j)
    {
      const size_t J = indexToCompare[j];

      const features::Regions *regionsJ = regions_perImage.at(J).get();
      const size_t regions_countJ = regionsJ->RegionCount();
      const typename MatcherT::ScalarT * tabJ = NULL;
      if(regions_countJ > 0)
        tabJ = reinterpret_cast<const typename MatcherT::ScalarT *>(regionsJ->DescriptorRawData());

      const size_t NNN__ = 2;
      std::vector<int> vec_nIndice10;
      std::vector<typename MatcherT::DistanceType> vec_fDistance10;

      //Find left->right
      matcher10.SearchNeighbours(tabJ, regions_countJ, &vec_nIndice10, &vec_fDistance10, NNN__);

      std::vector<IndMatch> vec_FilteredMatches;
      std::vector<int> vec_NNRatioIndexes;
      NNdistanceRatio( vec_fDistance10.begin(), // distance start
        vec_fDistance10.end(),  // distance end
        NNN__, // Number of neighbor in iterator sequence (minimum required 2)
        vec_NNRatioIndexes, // output (index that respect Lowe Ratio)
        fDistRatio);

      for (size_t k=0; k < vec_NNRatioIndexes.size(); ++k)
      {
        const size_t index = vec_NNRatioIndexes[k];
        vec_FilteredMatches.push_back(
          IndMatch(vec_nIndice10[index*NNN__], index) );
      }

      // Remove duplicates
      IndMatch::getDeduplicated(vec_FilteredMatches);

      // Remove matches that have the same (X,Y) coordinates
      const std::vector<PointFeature> pointFeaturesJ = regionsJ->GetRegionsPositions();
      IndMatchDecorator<float> matchDeduplicator(vec_FilteredMatches, pointFeaturesI, pointFeaturesJ);
      matchDeduplicator.getDeduplicated(vec_FilteredMatches);

      {
        ++my_progress_bar;
        if (!vec_FilteredMatches.empty())
        {
          map_PutativesMatches.insert( make_pair( make_pair(I,J), std::move(vec_FilteredMatches) ));
        }
      }
    }
  }
}

void Matcher_Regions_AllInMemory::Match(
  const std::vector<std::string> & vec_fileNames, // input filenames,
  const Pair_Set & pairs,
  PairWiseMatches & map_PutativesMatches)const // the pairwise photometric corresponding points
{
  if (regions_perImage.size() < 2)
  {
    return; // No sufficient images to compare (nothing to do)
  }
  else
  {
    // Build the required abstract Matchers from the regions Types
    const features::Regions *regions = regions_perImage.begin()->second.get();

    // Handle invalid request
    if (regions->IsScalar() && _eMatcherType == BRUTE_FORCE_HAMMING)
      return ;
    if (regions->IsBinary() && _eMatcherType != BRUTE_FORCE_HAMMING)
      return ;

    // Switch regions type ID, matcher & Metric: call the good MatcherT
    if (regions->IsScalar())
    {
      if(regions->Type_id() == typeid(unsigned char).name())
      {
        // Build on the fly unsigned char based Matcher
        switch(_eMatcherType)
        {
          case BRUTE_FORCE_L2:
          {
            typedef L2_Vectorized<unsigned char> MetricT;
            typedef ArrayMatcherBruteForce<unsigned char, MetricT> MatcherT;
            /// Match the distRatio to the used metric
            Template_Matcher<MatcherT>(vec_fileNames, pairs, map_PutativesMatches,
              regions_perImage, Square(fDistRatio), _eMatcherType);
          }
          break;
          case ANN_L2:
          {
            typedef flann::L2<unsigned char> MetricT;
            typedef ArrayMatcher_Kdtree_Flann<unsigned char, MetricT> MatcherT;
            /// Match the distRatio to the used metric
            Template_Matcher<MatcherT>(vec_fileNames, pairs, map_PutativesMatches,
              regions_perImage, Square(fDistRatio), _eMatcherType);
          }
          break;
          default:
            std::cerr << "Using unknown matcher type" << std::endl;
        }
      }
      else
      if(regions->Type_id() == typeid(float).name())
      {
        // Build on the fly float based Matcher
        switch(_eMatcherType)
        {
          case BRUTE_FORCE_L2:
          {
            typedef L2_Vectorized<float> MetricT;
            typedef ArrayMatcherBruteForce<float, MetricT> MatcherT;
            Template_Matcher<MatcherT>(vec_fileNames, pairs, map_PutativesMatches,
              regions_perImage, Square(fDistRatio), _eMatcherType);
          }
          break;
          case ANN_L2:
          {
            typedef flann::L2<float> MetricT;
            typedef ArrayMatcher_Kdtree_Flann<float, MetricT> MatcherT;
            /// Match the distRatio to the used metric
            Template_Matcher<MatcherT>(vec_fileNames, pairs, map_PutativesMatches,
              regions_perImage, Square(fDistRatio), _eMatcherType);
          }
          break;
          default:
            std::cerr << "Using unknown matcher type" << std::endl;
        }
      }
      else
      if(regions->Type_id() == typeid(double).name())
      {
        // Build on the fly double based Matcher
        switch(_eMatcherType)
        {
          case BRUTE_FORCE_L2:
          {
            typedef L2_Vectorized<double> MetricT;
            typedef ArrayMatcherBruteForce<double, MetricT> MatcherT;
            /// Match the distRatio to the used metric
            Template_Matcher<MatcherT>(vec_fileNames, pairs, map_PutativesMatches,
              regions_perImage, Square(fDistRatio), _eMatcherType);
          }
          break;
          case ANN_L2:
          {
            typedef flann::L2<double> MetricT;
            typedef ArrayMatcher_Kdtree_Flann<double, MetricT> MatcherT;
            /// Match the distRatio to the used metric
            Template_Matcher<MatcherT>(vec_fileNames, pairs, map_PutativesMatches,
              regions_perImage, Square(fDistRatio), _eMatcherType);
          }
          break;
          default:
            std::cerr << "Using unknown matcher type" << std::endl;
        }
      }
    }
    else
    if (regions->IsBinary() && regions->Type_id() == typeid(unsigned char).name())
    {
      switch(_eMatcherType)
      {
        case BRUTE_FORCE_HAMMING:
        {
          typedef Hamming<unsigned char> Metric;
          typedef ArrayMatcherBruteForce<unsigned char, Metric> MatcherT;
          Template_Matcher<MatcherT>(vec_fileNames, pairs, map_PutativesMatches,
           regions_perImage, fDistRatio, _eMatcherType);
        }
        break;
        default:
            std::cerr << "Using unknown matcher type" << std::endl;
      }
    }
  }
}

}; // namespace openMVG
