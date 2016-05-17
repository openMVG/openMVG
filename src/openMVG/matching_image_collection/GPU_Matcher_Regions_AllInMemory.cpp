
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/matching_image_collection/Cascade_Hashing_Matcher_Regions_AllInMemory.hpp"
#include "openMVG/matching/matcher_cascade_hashing.hpp"
#include "openMVG/matching/indMatchDecoratorXY.hpp"
#include "openMVG/matching/matching_filters.hpp"

#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "third_party/progress/progress.hpp"

namespace openMVG {
namespace matching_image_collection {

using namespace openMVG::matching;
using namespace openMVG::features;

Cascade_Hashing_Matcher_Regions_AllInMemory
::Cascade_Hashing_Matcher_Regions_AllInMemory
(
  float distRatio
):Matcher(), f_dist_ratio_(distRatio)
{
}

namespace impl
{
template <typename ScalarT>
void Match
(
  const sfm::SfM_Data & sfm_data,
  const sfm::Regions_Provider & regions_provider,
  const Pair_Set & pairs,
  float fDistRatio,
  PairWiseMatches & map_PutativesMatches // the pairwise photometric corresponding points
)
{
  C_Progress_display my_progress_bar( pairs.size() );

  // Collect used view indexes
  std::set<IndexT> used_index;
  // Sort pairs according the first index to minimize later memory swapping
  typedef std::map<IndexT, std::vector<IndexT> > Map_vectorT;
  Map_vectorT map_Pairs;
  for (Pair_Set::const_iterator iter = pairs.begin(); iter != pairs.end(); ++iter)
  {
    map_Pairs[iter->first].push_back(iter->second);
    used_index.insert(iter->first);
    used_index.insert(iter->second);
  }

  typedef Eigen::Matrix<ScalarT, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> BaseMat;

  // Init the cascade hasher
  CascadeHasher cascade_hasher;
  if (!used_index.empty())
  {
    const IndexT I = *used_index.begin();
    const features::Regions &regionsI = *regions_provider.regions_per_view.at(I).get();
    const size_t dimension = regionsI.DescriptorLength();
    cascade_hasher.Init(dimension);
  }

  std::map<IndexT, HashedDescriptions> hashed_base_;

  // Compute the zero mean descriptor that will be used for hashing (one for all the image regions)
  Eigen::VectorXf zero_mean_descriptor;
  {
    Eigen::MatrixXf matForZeroMean;
    for (int i =0; i < used_index.size(); ++i)
    {
      std::set<IndexT>::const_iterator iter = used_index.begin();
      std::advance(iter, i);
      const IndexT I = *iter;
      const features::Regions &regionsI = *regions_provider.regions_per_view.at(I).get();
      const ScalarT * tabI =
        reinterpret_cast<const ScalarT*>(regionsI.DescriptorRawData());
      const size_t dimension = regionsI.DescriptorLength();
      if (i==0)
      {
        matForZeroMean.resize(used_index.size(), dimension);
        matForZeroMean.fill(0.0f);
      }
      if (regionsI.RegionCount() > 0)
      {
        Eigen::Map<BaseMat> mat_I( (ScalarT*)tabI, regionsI.RegionCount(), dimension);
        matForZeroMean.row(i) = CascadeHasher::GetZeroMeanDescriptor(mat_I);
      }
    }
    zero_mean_descriptor = CascadeHasher::GetZeroMeanDescriptor(matForZeroMean);
  }

  // Index the input regions
#ifdef OPENMVG_USE_OPENMP
  #pragma omp parallel for schedule(dynamic)
#endif
  for (int i =0; i < used_index.size(); ++i)
  {
    std::set<IndexT>::const_iterator iter = used_index.begin();
    std::advance(iter, i);
    const IndexT I = *iter;
    const features::Regions &regionsI = *regions_provider.regions_per_view.at(I).get();
    const ScalarT * tabI =
      reinterpret_cast<const ScalarT*>(regionsI.DescriptorRawData());
    const size_t dimension = regionsI.DescriptorLength();

    Eigen::Map<BaseMat> mat_I( (ScalarT*)tabI, regionsI.RegionCount(), dimension);
    HashedDescriptions hashed_description = cascade_hasher.CreateHashedDescriptions(mat_I,
      zero_mean_descriptor);
#ifdef OPENMVG_USE_OPENMP
    #pragma omp critical
#endif
    {
      hashed_base_[I] = std::move(hashed_description);
    }
  }

  // Perform matching between all the pairs
  for (Map_vectorT::const_iterator iter = map_Pairs.begin();
    iter != map_Pairs.end(); ++iter)
  {
    const IndexT I = iter->first;
    const std::vector<IndexT> & indexToCompare = iter->second;

    const features::Regions &regionsI = *regions_provider.regions_per_view.at(I).get();
    if (regionsI.RegionCount() == 0)
    {
      my_progress_bar += indexToCompare.size();
      continue;
    }

    const std::vector<features::PointFeature> pointFeaturesI = regionsI.GetRegionsPositions();
    const ScalarT * tabI =
      reinterpret_cast<const ScalarT*>(regionsI.DescriptorRawData());
    const size_t dimension = regionsI.DescriptorLength();
    Eigen::Map<BaseMat> mat_I( (ScalarT*)tabI, regionsI.RegionCount(), dimension);

#ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for schedule(dynamic)
#endif
    for (int j = 0; j < (int)indexToCompare.size(); ++j)
    {
      size_t J = indexToCompare[j];
      const features::Regions &regionsJ = *regions_provider.regions_per_view.at(J).get();

      if (regions_provider.regions_per_view.count(J) == 0
          || regionsI.Type_id() != regionsJ.Type_id())
      {
#ifdef OPENMVG_USE_OPENMP
        #pragma omp critical
#endif
        ++my_progress_bar;
        continue;
      }

      // Matrix representation of the query input data;
      const ScalarT * tabJ = reinterpret_cast<const ScalarT*>(regionsJ.DescriptorRawData());
      Eigen::Map<BaseMat> mat_J( (ScalarT*)tabJ, regionsJ.RegionCount(), dimension);

      IndMatches pvec_indices;
      typedef typename Accumulator<ScalarT>::Type ResultType;
      std::vector<ResultType> pvec_distances;
      pvec_distances.reserve(regionsJ.RegionCount() * 2);
      pvec_indices.reserve(regionsJ.RegionCount() * 2);

      // Match the query descriptors to the database
      cascade_hasher.Match_HashedDescriptions<BaseMat, ResultType>(
        hashed_base_[J], mat_J,
        hashed_base_[I], mat_I,
        &pvec_indices, &pvec_distances);

      std::vector<int> vec_nn_ratio_idx;
      // Filter the matches using a distance ratio test:
      //   The probability that a match is correct is determined by taking
      //   the ratio of distance from the closest neighbor to the distance
      //   of the second closest.
      matching::NNdistanceRatio(
        pvec_distances.begin(), // distance start
        pvec_distances.end(),   // distance end
        2, // Number of neighbor in iterator sequence (minimum required 2)
        vec_nn_ratio_idx, // output (indices that respect the distance Ratio)
        Square(fDistRatio));

      matching::IndMatches vec_putative_matches;
      vec_putative_matches.reserve(vec_nn_ratio_idx.size());
      for (size_t k=0; k < vec_nn_ratio_idx.size(); ++k)
      {
        const size_t index = vec_nn_ratio_idx[k];
        vec_putative_matches.emplace_back(pvec_indices[index*2].j_, pvec_indices[index*2].i_);
      }

      // Remove duplicates
      matching::IndMatch::getDeduplicated(vec_putative_matches);

      // Remove matches that have the same (X,Y) coordinates
      const std::vector<features::PointFeature> pointFeaturesJ = regionsJ.GetRegionsPositions();
      matching::IndMatchDecorator<float> matchDeduplicator(vec_putative_matches,
        pointFeaturesI, pointFeaturesJ);
      matchDeduplicator.getDeduplicated(vec_putative_matches);

#ifdef OPENMVG_USE_OPENMP
#pragma omp critical
#endif
      {
        ++my_progress_bar;
        if (!vec_putative_matches.empty())
        {
          map_PutativesMatches.insert( make_pair( make_pair(I,J), std::move(vec_putative_matches) ));
        }
      }
    }
  }
}
} // namespace impl

void Cascade_Hashing_Matcher_Regions_AllInMemory::Match
(
  const sfm::SfM_Data & sfm_data,
  const std::shared_ptr<sfm::Regions_Provider> & regions_provider,
  const Pair_Set & pairs,
  PairWiseMatches & map_PutativesMatches // the pairwise photometric corresponding points
)const
{
#ifdef OPENMVG_USE_OPENMP
  std::cout << "Using the OPENMP thread interface" << std::endl;
#endif

  if (regions_provider->regions_per_view.empty())
    return;

  const features::Regions &regions = *regions_provider->regions_per_view.begin()->second.get();

  if (regions.IsBinary())
    return;

  if(regions.Type_id() == typeid(unsigned char).name())
  {
    impl::Match<unsigned char>(
      sfm_data,
      *regions_provider.get(),
      pairs,
      f_dist_ratio_,
      map_PutativesMatches);
  }
  else
  if(regions.Type_id() == typeid(float).name())
  {
    impl::Match<float>(
      sfm_data,
      *regions_provider.get(),
      pairs,
      f_dist_ratio_,
      map_PutativesMatches);
  }
  else
  {
    std::cerr << "Matcher not implemented for this region type" << std::endl;
  }
}

} // namespace openMVG
} // namespace matching_image_collection
