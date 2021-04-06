// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2020 Pierre MOULON, Romain JANVIER

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MATCHING_IMAGE_COLLECTION_VLAD_HPP
#define OPENMVG_MATCHING_IMAGE_COLLECTION_VLAD_HPP

#include "openMVG/matching_image_collection/VladBase.hpp"

#include "openMVG/clustering/kmeans.hpp"
#include "openMVG/features/descriptor.hpp"
#include "openMVG/matching/regions_matcher.hpp"

namespace openMVG {

template<typename RegionTypeT>
class VLAD : public VLADBase
{
  public:

  DescriptorVector RegionsToCodebook(
    const std::vector<IndexT>& view_ids,
    std::shared_ptr<sfm::Regions_Provider> learning_regions_provider
  ) override
  {
    using ScalarT = typename RegionTypeT::DescriptorT::bin_type;
    using ConstMatrixRef =
      Eigen::Map<const Eigen::Matrix<ScalarT, Eigen::Dynamic, Eigen::Dynamic,
                                     Eigen::RowMajor>>;

    DescriptorVector descriptor_array;

    const size_t base_descriptor_length =
        learning_regions_provider->getRegionsType()->DescriptorLength();

    for (const auto &view_id : view_ids) {
      const auto &regions = learning_regions_provider->get(view_id);
      RegionTypeT *cast_centroid_regions =
        dynamic_cast<RegionTypeT *>(regions.get());

      const ScalarT *tab =
          reinterpret_cast<const ScalarT *>(cast_centroid_regions->DescriptorRawData());
      ConstMatrixRef descriptors(tab, cast_centroid_regions->RegionCount(),
                                base_descriptor_length);
      for (int region_id = 0; region_id < cast_centroid_regions->RegionCount(); ++region_id) {
        descriptor_array.emplace_back(
            descriptors.row(region_id).template cast<typename DescriptorType::Scalar>());
      }
    }
    return descriptor_array;
  }

  DescriptorVector BuildCodebook(
    const DescriptorVector& descriptor_array,
    const int codebook_size = 128,
    const int max_nb_iteration = 25) override
  {
    DescriptorVector codebook;
    std::vector<uint32_t> vec_ids;
    const clustering::KMeansInitType k_mean_init_type =
        clustering::KMeansInitType::
            KMEANS_INIT_RANDOM;  // Kmean PP is way too slow vs. Random init (even
                                // if PP is better for corner cases)

    const bool do_progress_display = max_nb_iteration != std::numeric_limits<uint32_t>::max();
    system::LoggerProgress progress;
    clustering::KMeans(
        descriptor_array,
        vec_ids,
        codebook,
        codebook_size,
        max_nb_iteration,
        k_mean_init_type,
        &progress);
    return codebook;
  }

  // Utility function
  void CodebookToRegions(
    std::unique_ptr<features::Regions>& centroid_regions,
    const DescriptorVector& codebook
  ) const override
  {
    RegionTypeT *cast_centroid_regions =
        dynamic_cast<RegionTypeT *>(centroid_regions.get());
    cast_centroid_regions->Descriptors().reserve(codebook.size());
    cast_centroid_regions->Features().resize(codebook.size());

    for (const auto &centroid : codebook) {
      cast_centroid_regions->Descriptors().emplace_back(
          centroid.template cast<typename RegionTypeT::DescriptorT::bin_type>());
    }
    // Note: We ignore filling feature data since we don't need them
  }

  VladMatrixType ComputeVLADEmbedding(
    const std::vector<IndexT>& view_ids,
    std::unique_ptr<features::Regions>& centroid_regions, // The codebook
    std::shared_ptr<sfm::Regions_Provider> embedding_regions_provider,
    const VLAD_NORMALIZATION vlad_normalization_type =
    VLAD_NORMALIZATION::RESIDUAL_NORMALIZATION_PWR_LAW)
    override
  {
    system::LoggerProgress progress;

    const size_t codebook_size = centroid_regions->RegionCount();
    const size_t base_descriptor_length =
        embedding_regions_provider->getRegionsType()->DescriptorLength();
    const size_t vlad_descriptor_length = base_descriptor_length * codebook_size;

    VladMatrixType mat_vlad_descriptors =
      VladMatrixType::Zero(vlad_descriptor_length, view_ids.size());
    // For each image (regions), compute its VLAD representation
    const RegionTypeT *cast_centroid_regions =
        dynamic_cast<RegionTypeT *>(centroid_regions.get());
    Vec vlad_desc(vlad_descriptor_length);
    progress.Restart(
        view_ids.size(), "- VLAD Embedding... -");
    for (const auto &view_id : view_ids) {
      vlad_desc.setZero();
      const auto &query_regions = embedding_regions_provider->get(view_id);

      // Retrieve indexes (match descriptors to centroids)
      matching::IndMatches centroid_to_descriptor_associations;
      matching::Match(matching::EMatcherType::BRUTE_FORCE_L2, *centroid_regions,
                      *query_regions, centroid_to_descriptor_associations);

      const RegionTypeT *cast_query_regions =
          dynamic_cast<RegionTypeT *>(query_regions.get());

      // Accumulation of residual to the centroid
      for (const auto centroid_id_and_descriptor_list :
          centroid_to_descriptor_associations) {
        const auto centroid_id = centroid_id_and_descriptor_list.i_;
        const auto descriptor_id = centroid_id_and_descriptor_list.j_;

        const auto residual =
            cast_query_regions->Descriptors()[descriptor_id].template cast<double>() -
            cast_centroid_regions->Descriptors()[centroid_id].template cast<double>();

        switch (vlad_normalization_type) {
          case VLAD_NORMALIZATION::RESIDUAL_NORMALIZATION_PWR_LAW: {
          const auto norm_residual = residual.normalized();
          vlad_desc.segment(centroid_id * base_descriptor_length,
                            base_descriptor_length) += residual;
        }
        break;
        default:
          vlad_desc.segment(centroid_id * base_descriptor_length,
                            base_descriptor_length) += residual;
        }
      }

      // per descriptor normalization
      for (int centroid_id = 0; centroid_id < codebook_size; ++centroid_id) {
        auto local_vlad = vlad_desc.segment(centroid_id * base_descriptor_length,
                                            base_descriptor_length);
        switch (vlad_normalization_type) {
          case VLAD_NORMALIZATION::INTRA_NORMALIZATION:
            local_vlad.normalize();
            break;
          case VLAD_NORMALIZATION::SIGNED_SQUARE_ROOTING:
            local_vlad =
                local_vlad.array().sign() * local_vlad.array().abs().sqrt();
            break;
          case VLAD_NORMALIZATION::RESIDUAL_NORMALIZATION_PWR_LAW:
            local_vlad =
                local_vlad.array().sign() * local_vlad.array().abs().pow(0.2);
            break;
        }
      }

      // if(max_feats > 0) TODO(RJ): center adaptation for All About VLAD

      // Global L2 normalization, it is used by all variants of VLAD
      vlad_desc.normalize();

      // Insert the vector into the matrix
      mat_vlad_descriptors.col(view_id) =
          vlad_desc.cast<VladMatrixType::Scalar>();
      ++progress;
    }
    return mat_vlad_descriptors;
  }
};


} // namespace openMVG

#endif // OPENMVG_MATCHING_IMAGE_COLLECTION_VLAD_HPP
