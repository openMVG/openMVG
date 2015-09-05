
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
#include "openMVG/matching/regions_matcher.hpp"

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
    const std::vector<size_t> & indexToCompare = iter->second;

    if (regions_provider->regions_per_view.count(I) == 0)
      continue;

    const features::Regions& regionsI = *regions_provider->regions_per_view.at(I).get();
    if (regionsI.RegionCount() == 0)
    {
      my_progress_bar += indexToCompare.size();
      continue;
    }
    
    matching::RegionsMatcher<MatcherT> matcher(regionsI);

    for (int j = 0; j < (int)indexToCompare.size(); ++j, ++my_progress_bar)
    {
      const size_t J = indexToCompare[j];

      if (regions_provider->regions_per_view.count(J) == 0)
        continue;

      const features::Regions& regionsJ = *regions_provider->regions_per_view.at(J).get();
      if(regionsJ.RegionCount() == 0)
        continue;
      
      std::vector<IndMatch> vec_matches;
      matcher.MatchRatioTest(vec_matches, regionsJ, fDistRatio);

      if (!vec_matches.empty())
      {
        map_PutativesMatches.insert( make_pair( make_pair(I,J), std::move(vec_matches) ));
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
