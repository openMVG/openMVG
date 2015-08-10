
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/matching_image_collection/Matcher_Regions_AllInMemory.hpp"

#include "openMVG/features/regions.hpp"
#include "openMVG/matching_image_collection/Matcher.hpp"
#include "openMVG/matching/matcher_brute_force.hpp"
#include "openMVG/matching/matcher_kdtree_flann.hpp"
#include "openMVG/matching/matcher_cascade_hashing.hpp"
#include "openMVG/matching/indMatchDecoratorXY.hpp"
#include "openMVG/matching/matching_filters.hpp"

#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "third_party/progress/progress.hpp"

namespace openMVG {
namespace matching_image_collection {

using namespace openMVG::matching;
using namespace openMVG::features;

Matcher_Regions_AllInMemory::Matcher_Regions_AllInMemory(
  float distRatio, EMatcherType eMatcherType)
  :Matcher(), fDistRatio(distRatio), _eMatcherType(eMatcherType)
{
}

/// Template matching
template <typename MatcherT>
void Template_Matcher(
  const Pair_Set & pairs,
  PairWiseMatches & map_PutativesMatches, // the pairwise photometric corresponding points
  // Data & parameters
  const std::shared_ptr<sfm::Regions_Provider> & regions_provider,
  // Distance ratio used to discard spurious correspondence
  const float fDistRatio,
  // Matcher Type
  EMatcherType _eMatcherType)
{
#ifdef OPENMVG_USE_OPENMP
  std::cout << "Using the OPENMP thread interface" << std::endl;
#endif
  bool b_multithreaded_pair_search = false;
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
    case CASCADE_HASHING_L2:
      std::cout << "Using CASCADE_HASHING_L2 matcher" << std::endl;
      b_multithreaded_pair_search = true;
      // -> set to true only here, since OpenMP instructions are not used in this matcher
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
    const std::vector<size_t> & indexToCompare = iter->second;

    if (regions_provider->regions_per_view.count(I) == 0)
      continue;

    const features::Regions *regionsI = regions_provider->regions_per_view.at(I).get();
    const size_t regions_countI = regionsI->RegionCount();
    if (regions_countI == 0)
    {
      my_progress_bar += indexToCompare.size();
      continue;
    }

    const std::vector<PointFeature> pointFeaturesI = regionsI->GetRegionsPositions();
    const typename MatcherT::ScalarT * tabI =
      reinterpret_cast<const typename MatcherT::ScalarT *>(regionsI->DescriptorRawData());

    MatcherT matcher10;
    if (!matcher10.Build(tabI, regions_countI, regionsI->DescriptorLength()))
    {
      my_progress_bar += indexToCompare.size();
      continue;
    }
#ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for schedule(dynamic) if(b_multithreaded_pair_search)
#endif
    for (int j = 0; j < (int)indexToCompare.size(); ++j)
    {
      const size_t J = indexToCompare[j];

      if (regions_provider->regions_per_view.count(J) == 0)
        continue;

      const features::Regions *regionsJ = regions_provider->regions_per_view.at(J).get();
      const size_t regions_countJ = regionsJ->RegionCount();
      if(regions_countJ == 0)
      {
        continue;
      }

      const typename MatcherT::ScalarT * tabJ =
        reinterpret_cast<const typename MatcherT::ScalarT *>(regionsJ->DescriptorRawData());

      const size_t NNN__ = 2;
      std::vector<int> vec_nIndice10;
      std::vector<typename MatcherT::DistanceType> vec_fDistance10;

      //Find left->right
      if (matcher10.SearchNeighbours(tabJ, regions_countJ, &vec_nIndice10, &vec_fDistance10, NNN__))
      {
        std::vector<IndMatch> vec_FilteredMatches;
        std::vector<int> vec_NNRatioIndexes;
        NNdistanceRatio( vec_fDistance10.begin(), // distance start
          vec_fDistance10.end(),  // distance end
          NNN__, // Number of neighbor in iterator sequence (minimum required 2)
          vec_NNRatioIndexes, // output (indices that respect Lowe Ratio)
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

#ifdef OPENMVG_USE_OPENMP
  #pragma omp critical
#endif
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
}

void Matcher_Regions_AllInMemory::Match(
  const sfm::SfM_Data & sfm_data,
  const std::shared_ptr<sfm::Regions_Provider> & regions_provider,
  const Pair_Set & pairs,
  PairWiseMatches & map_PutativesMatches)const // the pairwise photometric corresponding points
{
  if (regions_provider->regions_per_view.size() < 2)
  {
    return; // No sufficient images to compare (nothing to do)
  }
  else
  {
    // Build the required abstract Matchers according of the regions descriptor Types
    const features::Regions *regions = regions_provider->regions_per_view.begin()->second.get();

    // Handle invalid request
    if (regions->IsScalar() && _eMatcherType == BRUTE_FORCE_HAMMING)
      return ;
    if (regions->IsBinary() && _eMatcherType != BRUTE_FORCE_HAMMING)
      return ;

    // Switch regions type ID, matcher & Metric: call the good MatcherT
    if (regions->IsScalar())
    {
      if (regions->Type_id() == typeid(unsigned char).name())
      {
        // Build on the fly unsigned char based Matcher
        switch (_eMatcherType)
        {
          case BRUTE_FORCE_L2:
          {
            typedef L2_Vectorized<unsigned char> MetricT;
            typedef ArrayMatcherBruteForce<unsigned char, MetricT> MatcherT;
            /// Match the distRatio to the used metric
            Template_Matcher<MatcherT>(pairs, map_PutativesMatches,
              regions_provider, Square(fDistRatio), _eMatcherType);
          }
          break;
          case ANN_L2:
          {
            typedef flann::L2<unsigned char> MetricT;
            typedef ArrayMatcher_Kdtree_Flann<unsigned char, MetricT> MatcherT;
            /// Match the distRatio to the used metric
            Template_Matcher<MatcherT>(pairs, map_PutativesMatches,
              regions_provider, Square(fDistRatio), _eMatcherType);
          }
          break;
          case CASCADE_HASHING_L2:
          {
            typedef L2_Vectorized<unsigned char> MetricT;
            typedef ArrayMatcherCascadeHashing<unsigned char, MetricT> MatcherT;
            /// Match the distRatio to the used metric
            Template_Matcher<MatcherT>(pairs, map_PutativesMatches,
              regions_provider, Square(fDistRatio), _eMatcherType);
          }
          break;
          default:
            std::cerr << "Using unknown matcher type" << std::endl;
        }
      }
      else
      if (regions->Type_id() == typeid(float).name())
      {
        // Build on the fly float based Matcher
        switch (_eMatcherType)
        {
          case BRUTE_FORCE_L2:
          {
            typedef L2_Vectorized<float> MetricT;
            typedef ArrayMatcherBruteForce<float, MetricT> MatcherT;
            Template_Matcher<MatcherT>(pairs, map_PutativesMatches,
              regions_provider, Square(fDistRatio), _eMatcherType);
          }
          break;
          case ANN_L2:
          {
            typedef flann::L2<float> MetricT;
            typedef ArrayMatcher_Kdtree_Flann<float, MetricT> MatcherT;
            /// Match the distRatio to the used metric
            Template_Matcher<MatcherT>(pairs, map_PutativesMatches,
              regions_provider, Square(fDistRatio), _eMatcherType);
          }
          break;
          case CASCADE_HASHING_L2:
          {
            typedef L2_Vectorized<float> MetricT;
            typedef ArrayMatcherCascadeHashing<float, MetricT> MatcherT;
            /// Match the distRatio to the used metric
            Template_Matcher<MatcherT>(pairs, map_PutativesMatches,
              regions_provider, Square(fDistRatio), _eMatcherType);
          }
          break;
          default:
            std::cerr << "Using unknown matcher type" << std::endl;
        }
      }
      else
      if (regions->Type_id() == typeid(double).name())
      {
        // Build on the fly double based Matcher
        switch (_eMatcherType)
        {
          case BRUTE_FORCE_L2:
          {
            typedef L2_Vectorized<double> MetricT;
            typedef ArrayMatcherBruteForce<double, MetricT> MatcherT;
            /// Match the distRatio to the used metric
            Template_Matcher<MatcherT>(pairs, map_PutativesMatches,
              regions_provider, Square(fDistRatio), _eMatcherType);
          }
          break;
          case ANN_L2:
          {
            typedef flann::L2<double> MetricT;
            typedef ArrayMatcher_Kdtree_Flann<double, MetricT> MatcherT;
            /// Match the distRatio to the used metric
            Template_Matcher<MatcherT>(pairs, map_PutativesMatches,
              regions_provider, Square(fDistRatio), _eMatcherType);
          }
          break;
          case CASCADE_HASHING_L2:
          {
            std::cerr << "Not yet implemented" << std::endl;
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
      switch (_eMatcherType)
      {
        case BRUTE_FORCE_HAMMING:
        {
          typedef Hamming<unsigned char> Metric;
          typedef ArrayMatcherBruteForce<unsigned char, Metric> MatcherT;
          Template_Matcher<MatcherT>(pairs, map_PutativesMatches,
           regions_provider, fDistRatio, _eMatcherType);
        }
        break;
        default:
            std::cerr << "Using unknown matcher type" << std::endl;
      }
    }
    else
    {
      std::cerr << "Please consider add this region type_id to Matcher_Regions_AllInMemory::Match(...)\n"
        << "typeid: " << regions->Type_id() << std::endl;
    }
  }
}

} // namespace openMVG
} // namespace matching_image_collection
