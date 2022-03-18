// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 JiaWang BIAN, Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_ROBUST_ESTIMATION_GMSFILTER_HPP
#define OPENMVG_ROBUST_ESTIMATION_GMSFILTER_HPP

#include "openMVG/matching/indMatch.hpp"
#include "openMVG/numeric/eigen_alias_definition.hpp"

#include <array>
#include <map>
#include <utility>

namespace openMVG {
namespace robust {

//  Implementation of the GMS Filter (Grid-based Motion Statistics filter) [1]
// [1] GMS: Grid-based Motion Statistics for Fast, Ultra-robust Feature Correspondence.
//     JiaWang Bian, Wen-Yan Lin, Yasuyuki Matsushita, Sai-Kit Yeung, Tan Dat Nguyen, Ming-Ming Cheng
//     CVPR, 2017

class GMSFilter
{
public:
  /// @brief GMS constructor (setup corresponding point positions & image sizes)
  /// @param [in] point_positions1    Image 1 point positions
  /// @param [in] image_size1         Image 1 size (w,h)
  /// @param [in] point_positions2    Image 2 point positions
  /// @param [in] image_size2         Image 2 size (w,h)
  /// @param [in] matches             Corresponding indexes for point_position
  /// array that depict a match.
  /// @param [in] threshold_factor    GMS threshold (Original paper uses t=6,
  ///  We advise 12 in order to limit the number of false positives)
  GMSFilter
  (
    const std::vector<Eigen::Vector2f> & point_positions1,
    const std::pair<int,int> & image_size1,
    const std::vector<Eigen::Vector2f> & point_positions2,
    const std::pair<int,int> & image_size2,
    const matching::IndMatches & matches,
    const int threshold_factor = 12
  );

  ///  @brief Get the inliers indexes thanks to an inlier boolean mask
  ///
  ///  @param [out] inlier_mask         The inlier/outlier classification as a boolean mask
  ///  @param [in] scale_invariance     Ask for rotation invariance (Use many grid rotations)
  ///  @param [in] rotation_invariance  Ask for scale invariance
  ///
  ///  @return The number of found correspondences (inliers)
  int GetInlierMask
  (
    std::vector<bool> & inlier_mask,
    const bool scale_invariance = false,
    const bool rotation_invariance = false
  );

private:

  ///  @brief Normalize point coordinates to the range (0 - 1)
  ///
  ///  @param [in] image_size   The image size (w, h)
  ///  @param [in,out] points   The points that are normalized inplace
  ///
  void NormalizePoints
  (
    const std::pair<int, int> & image_size,
    std::vector<Eigen::Vector2f> & points
  ) const;

  ///  @brief Compute the grid index on the Left image for the given point and the chosen grid type
  ///
  ///  @param [in] pt     The normalized point position
  ///  @param [in] type   The grid type (1:center, 2: East, 3:South, 4:West)
  ///
  ///  @return The left grid index
  int GetGridIndexLeft
  (
    const Eigen::Vector2f & pt,
    int type
  ) const;

  ///  @brief Compute the grid index on the Right image for the given point and the chosen grid type
  ///
  ///  @param [in] pt     The normalized point position
  ///  @param [in] type   The subgrid type (1:center, 2:East, 3:South, 4:SEst)
  ///
  ///  @return The right grid index
  int GetGridIndexRight
  (
    const Eigen::Vector2f & pt
  ) const;

  ///  @brief Accumulate motion statistics for the match_pairs and given grid_type
  ///
  ///  @param [in] grid_type The grid type (1:center, 2:East, 3:South, 4:SEst)
  ///
  void AssignMatchPairs(int grid_type);

  ///  @brief Verify Cell Pairs
  ///  Threshold the aggregate statistics for the neighborhood grids.
  ///
  ///  @param [in] rotation_type The rotation type [0, 7]
  ///
  void VerifyCellPairs(int rotation_type);

  ///  @brief Get the 9 neighbors grid index
  ///
  ///  @param [in] idx       The grid index
  ///  @param [in] grid_size The grid size
  ///
  std::array<int, 9> GetNB9
  (
    const int idx,
    const std::pair<int, int> & grid_size
  ) const;

  ///  @brief Run the GMS filter for a given rotation_type and scale
  ///
  ///  @param [in] rotation_type  The rotation_type index [0, 7]
  ///  @param [in] scale          The scale index [0,4]
  ///
  int Run
  (
    int rotation_type,
    int scale
  );

  ///  @brief Compute the grid neighbor indexes
  ///
  ///  @param [in] neighborhood_grid_indexes   The image size (w, h)
  ///  @param [in,out] grid_size               The grid size
  ///
  void InitalizeGridNeighborhood
  (
    std::vector<std::array<int,9>> & neighborhood_grid_indexes,
    const std::pair<int, int> & grid_size
  ) const;

  ///  @brief Compute the Right grid  indexes according a scale index
  ///
  ///  @param [in,out] grid_size  The scale index [0, 4]
  ///
  void InitalizeRightGridNeighborhoodWithScale(int scale);

  // --
  // Data
  // --

  const uint8_t kGridSize = 20;

  // Normalized Points
  std::vector<Eigen::Vector2f> point_positions1_, point_positions2_;

  // Matches
  const matching::IndMatches matches_indexes_;

  // Threshold used to classify inlier / outlier motion statistics
  const int threshold_factor_;

  // Grid Size
  std::pair<int,int> left_grid_size_, right_grid_size_;
  // Left Grid count
  const int left_grid_count_;
  // Right Grid count
  int right_grid_count_;

  // x     : left grid idx
  // y     : right grid idx
  // value : how many matches from idx_left to idx_right
  using Mati = Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>;
  Mati motion_statistics_;

  // Number of point per cell
  std::vector<int> nb_points_per_left_cell_;

  // Index  : grid_idx_left
  // Value  : grid_idx_right
  std::vector<int> cell_pairs_;

  // Every Matches has a cell-pair
  // first  : grid_idx_left
  // second : grid_idx_right
  std::vector<std::pair<int, int>> match_cell_pairs_;

  // Inlier Mask (GMS kept correspondences)
  std::vector<bool> inlier_mask_;

  // Neighborhood grid indexes
  std::vector<std::array<int,9>>
    left_neighbor_grid_indexes_, right_neighbor_grid_indexes_;
};

} // namespace robust
} // namespace openMVG
#endif // OPENMVG_ROBUST_ESTIMATION_GMSFILTER_HPP
