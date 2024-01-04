// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2020 Pierre MOULON, Romain JANVIER

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MATCHING_IMAGE_COLLECTION_VLADBASE_HPP
#define OPENMVG_MATCHING_IMAGE_COLLECTION_VLADBASE_HPP

#include "openMVG/sfm/pipelines/sfm_regions_provider.hpp"

namespace openMVG {

// Vector of Locally Aggregated Descriptors (VLAD) possible vector normalization
enum class VLAD_NORMALIZATION : uint8_t {
  SIGNED_SQUARE_ROOTING =
      0,  // [1] Element wise power law normalization with alpha = 0.5 (SSR)
  INTRA_NORMALIZATION,  // [2] cluster wise L2 normalization, the so called
                        // "intra normalization"
  RESIDUAL_NORMALIZATION_PWR_LAW  // [3] Per residual L2 normalization (RN)
                                  // followed by and element wise power law
                                  // normalization with alpha = 0.2
};

// An implementation of Vector of Locally Aggregated Descriptors (VLAD) encoding
//
// [1]
// "Aggregating local descriptors into compact codes". H. Jegou, F.
// Perronnin, M. Douze, J. Sanchez, and P. Perez., C. Schmid. PAMI, 2012.
// [2]
// "All About VLAD". R. Arandjelovic and A. Zisserman. CVPR 2013.
// [3]
// "Revisiting the VLAD image representation". J. Delhumeau, P.H. Gosselin, H.
// Jégou, P. and Pérez. ACM Multimedia 2013.

class VLADBase
{
public:
  using KmeanInternalType = float;
  using VladInternalType = float;
  using DescriptorType = Eigen::Matrix<KmeanInternalType, Eigen::Dynamic, 1>;
  using DescriptorVector = std::vector<DescriptorType>;

  // VLADs
  using VladMatrixType =
        Eigen::Matrix<VladInternalType, Eigen::Dynamic, Eigen::Dynamic>;

  // IO
  // Convert Regions to a contiguous feature vector
  virtual DescriptorVector RegionsToCodebook(
    const std::vector<IndexT>& view_ids,
    std::shared_ptr<sfm::Regions_Provider> learning_regions_provider
  ) = 0;

  // IO
  // Convert a codebook to feature regions
  virtual void CodebookToRegions(
    std::unique_ptr<features::Regions>& centroid_regions,
    const DescriptorVector& codebook
  ) const = 0;


  // Build a codebook from a selection of descriptors
  virtual DescriptorVector BuildCodebook(
    const DescriptorVector& descriptor_array,
    const int codebook_size = 128,
    const int max_nb_iteration = 25) = 0;

  // Compute the VLAD representation of each "image" given the codebook
  // and its associated image descriptors
  virtual VladMatrixType ComputeVLADEmbedding(
    const std::vector<IndexT>& view_ids,
    std::unique_ptr<features::Regions>& centroid_regions, // The codebook
    std::shared_ptr<sfm::Regions_Provider> embedding_regions_provider,
    const VLAD_NORMALIZATION vlad_normalization_type =
      VLAD_NORMALIZATION::RESIDUAL_NORMALIZATION_PWR_LAW) = 0;
};

} // namespace openMVG

#endif // OPENMVG_MATCHING_IMAGE_COLLECTION_VLADBASE_HPP
