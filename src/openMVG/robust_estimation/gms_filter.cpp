// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 JiaWang BIAN, Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <openMVG/robust_estimation/gms_filter.hpp>

namespace openMVG {
namespace robust {

// Store the 8 possible 3 x 3 grid rotation pattern
static const std::array<std::array<uint8_t,9>,8> kRotationPatterns =
{
  1,2,3,
  4,5,6,
  7,8,9,

  4,1,2,
  7,5,3,
  8,9,6,

  7,4,1,
  8,5,2,
  9,6,3,

  8,7,4,
  9,5,1,
  6,3,2,

  9,8,7,
  6,5,4,
  3,2,1,

  6,9,8,
  3,5,7,
  2,1,4,

  3,6,9,
  2,5,8,
  1,4,7,

  2,3,6,
  1,5,9,
  4,7,8
};

// 5 level scales
static const std::array<double, 5> kScaleRatios = {
  { 1.0, 1.0 / 2.0, 1.0 / sqrt(2.0), sqrt(2.0), 2.0 }
};

GMSFilter::GMSFilter
(
  const std::vector<Eigen::Vector2f> & point_positions1,
  const std::pair<int,int> & image_size1,
  const std::vector<Eigen::Vector2f> & point_positions2,
  const std::pair<int,int> & image_size2,
  const matching::IndMatches & matches,
  const int threshold_factor
):
  point_positions1_(point_positions1),
  point_positions2_(point_positions2),
  matches_indexes_(matches),
  // Grid initialization
  left_grid_size_(kGridSize, kGridSize),
  left_grid_count_(left_grid_size_.first * left_grid_size_.second),
  threshold_factor_(threshold_factor)
{
  // Normalize the points position to [{0,1};{0,1}]
  NormalizePoints(image_size1, point_positions1_);
  NormalizePoints(image_size2, point_positions2_);

  // Initialize the right grid neighbor indexes
  left_neighbor_grid_indexes_.resize(left_grid_count_);
  InitalizeGridNeighborhood(left_neighbor_grid_indexes_, left_grid_size_);
};

void GMSFilter::NormalizePoints
(
  const std::pair<int, int> & image_size,
  std::vector<Eigen::Vector2f> & points
) const
{
  for (auto & kp : points)
  {
    kp.x() /= static_cast<float>(image_size.first);
    kp.y() /= static_cast<float>(image_size.second);
  }
}

void GMSFilter::InitalizeGridNeighborhood
(
  std::vector<std::array<int,9>> & neighborhood_grid_indexes,
  const std::pair<int, int> & grid_size
) const
{
  for (int i = 0; i < static_cast<int>(neighborhood_grid_indexes.size()); ++i)
  {
    neighborhood_grid_indexes[i] = GetNB9(i, grid_size);
  }
}

void GMSFilter::InitalizeRightGridNeighborhoodWithScale
(
  int scale
)
{
  // Set Scale
  right_grid_size_.first  = left_grid_size_.first  * kScaleRatios[scale];
  right_grid_size_.second = left_grid_size_.second * kScaleRatios[scale];
  right_grid_count_ = right_grid_size_.first * right_grid_size_.second;
  right_neighbor_grid_indexes_.resize(right_grid_count_);
  InitalizeGridNeighborhood(right_neighbor_grid_indexes_, right_grid_size_);
}

int GMSFilter::GetGridIndexLeft
(
  const Eigen::Vector2f & pt,
  int type
) const
{
  // Moving the point according the asked direction
  // (1:center, 2:East, 3:South, 4:SEst)
  // C -> E
  // | \
  // S   SE
  const double x_offset = (type == 2 || type == 4) ? 0.5 : 0;
  const double y_offset = (type == 3 || type == 4) ? 0.5 : 0;

  const int x = std::floor(pt.x() * left_grid_size_.first  + x_offset);
  const int y = std::floor(pt.y() * left_grid_size_.second + y_offset);

  // Be sure that the point still belong to the grid bounds
  if (x >= left_grid_size_.first ||
      y >= left_grid_size_.second ||
      x < 0 ||
      y < 0)
  {
    return -1;
  }

  return x + y * left_grid_size_.first;
}

int GMSFilter::GetGridIndexRight
(
  const Eigen::Vector2f & pt
) const
{
  const int x = std::floor(pt.x() * right_grid_size_.first);
  const int y = std::floor(pt.y() * right_grid_size_.second);

  // Be sure that the point is not out of bounds
  if (x >= right_grid_size_.first ||
      y >= right_grid_size_.second ||
      x < 0 ||
      y < 0)
  {
    return -1;
  }

  return x + y * right_grid_size_.first;
}

std::array<int, 9> GMSFilter::GetNB9
(
  const int idx,
  const std::pair<int, int> & grid_size
) const
{
  std::array<int, 9> NB9;

  const int idx_x = idx % grid_size.first;
  const int idx_y = idx / grid_size.first;

  for (const int yi : {-1, 0, 1})
  {
    const int idx_yy = idx_y + yi;
    for (const int xi : {-1, 0, 1})
    {
      const int idx_xx = idx_x + xi;
      const int index_NB9 = 4 + xi + yi * 3;
      if (idx_xx < 0 || idx_xx >= grid_size.first ||
          idx_yy < 0 || idx_yy >= grid_size.second)
        NB9[index_NB9] = -1;
      else
        NB9[index_NB9] = idx_xx + idx_yy * grid_size.first;
    }
  }
  return NB9;
}

int GMSFilter::GetInlierMask
(
  std::vector<bool> & inlier_mask,
  const bool scale_invariance,
  const bool rotation_invariance
)
{
  int max_inlier = 0;

  const std::vector<int> scales_idx = scale_invariance ?
    std::vector<int>({0, 1, 2, 3, 4}) :
    std::vector<int>({0});

  const std::vector<int> rotations_idx = rotation_invariance ?
    std::vector<int>({1, 2, 3, 4, 5, 6, 7, 8}) :
    std::vector<int>({1});

  for (const auto & scale_it : scales_idx)
  {
    InitalizeRightGridNeighborhoodWithScale(scale_it);
    for (const auto & rotation_it : rotations_idx)
    {
      const int num_inlier = Run(rotation_it, scale_it);
      if (num_inlier > max_inlier)
      {
        inlier_mask = inlier_mask_;
        max_inlier = num_inlier;
      }
    }
  }
  return max_inlier;
}

void GMSFilter::AssignMatchPairs(int grid_type)
{
  for (size_t i = 0; i < matches_indexes_.size(); ++i)
  {
    const auto & lp = point_positions1_[matches_indexes_[i].i_];
    const auto & rp = point_positions2_[matches_indexes_[i].j_];

    const int lgidx = match_cell_pairs_[i].first = GetGridIndexLeft(lp, grid_type);
    if (lgidx == -1) continue;
    const int rgidx = match_cell_pairs_[i].second = GetGridIndexRight(rp);
    if (rgidx == -1) continue;

    ++motion_statistics_(lgidx,rgidx);
    ++nb_points_per_left_cell_[lgidx];
  }
}

void GMSFilter::VerifyCellPairs(int rotation_type)
{
  for (int i = 0; i < left_grid_count_; ++i)
  {
    if (motion_statistics_.row(i).sum() == 0)
    {
      cell_pairs_[i] = -1;
      continue;
    }

    int max_number (0);
    for (int j = 0; j < right_grid_count_; j++)
    {
      if (motion_statistics_(i,j) > max_number)
      {
        cell_pairs_[i] = j;
        max_number     = motion_statistics_(i,j);
      }
    }

    const int idx_grid_rt = cell_pairs_[i];

    const auto & NB9_lt = left_neighbor_grid_indexes_[i];
    const auto & NB9_rt = right_neighbor_grid_indexes_[idx_grid_rt];

    int score (0);
    double thresh (0.0);
    int numpair (0);

    // Aggregate statistics for the neighborhood grids
    for (const int j : {0, 1, 2, 3, 4, 5, 6, 7, 8})
    {
      const auto & rotation_permutations = kRotationPatterns[rotation_type - 1];
      const int ll = NB9_lt[j];
      const int rr = NB9_rt[rotation_permutations[j] - 1];
      if (ll == -1 || rr == -1)
        continue;

      score += motion_statistics_(ll, rr);
      thresh += nb_points_per_left_cell_[ll];
      ++numpair;
    }
    if (numpair != 0)
    {
      thresh = threshold_factor_ * std::sqrt(thresh / static_cast<double>(numpair));
    }

    // Discard the match if it have a lower threshold or not computed statistics
    if (score < thresh || numpair == 0)
    {
      cell_pairs_[i] = -2;
    }
  }
}

int GMSFilter::Run
(
  int rotation_type,
  int scale
)
{
  const auto number_of_matches = matches_indexes_.size();

  inlier_mask_.assign(number_of_matches, false);

  // Initialize motion statistics
  match_cell_pairs_.assign(number_of_matches, {0, 0});
  motion_statistics_ = Mati::Zero(left_grid_count_, right_grid_count_);

  for (const int grid_type_it : {1, 2, 3, 4})
  {
    // Initialize arrays
    motion_statistics_.fill(0);
    cell_pairs_.assign(left_grid_count_, -1);
    nb_points_per_left_cell_.assign(left_grid_count_, 0);

    // Compute and verify motion statistics
    AssignMatchPairs(grid_type_it);
    VerifyCellPairs(rotation_type);

    // Mark inliers
    for (size_t i = 0; i < number_of_matches; ++i)
    {
      if (match_cell_pairs_[i].first  >= 0 &&
          match_cell_pairs_[i].second >= 0 &&
          cell_pairs_[match_cell_pairs_[i].first] == match_cell_pairs_[i].second)
      {
        inlier_mask_[i] = true;
      }
    }
  }
  // Return the number of inliers
  return std::count(inlier_mask_.cbegin(), inlier_mask_.cend(), true);
}

} // namespace robust
} // namespace openMVG
