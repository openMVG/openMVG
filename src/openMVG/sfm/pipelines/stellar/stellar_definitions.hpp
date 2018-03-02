// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2018 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_STELLAR_STELLAR_DEFINITIONS_HPP
#define OPENMVG_SFM_STELLAR_STELLAR_DEFINITIONS_HPP

#include "openMVG/sfm/pipelines/stellar/relative_scale.hpp"
#include "openMVG/numeric/l1_solver_admm.hpp"

#include <Eigen/Sparse>

namespace openMVG{
namespace sfm{

/// Mode that can be used to solve to rescale the relative translation of a stellar pod
enum class Stellar_Translation_Averaging_Solver_Type
{
  SCALING_SOLVER_L1,
  SCALING_SOLVER_L2,
  SCALING_SOLVER_L2_FULL
};

/// Check if a stellar pod defined by some relative scale is one CC or not
inline bool Relative_scales_are_one_cc
(
  const std::vector<Relative_Scale> & relative_scales
)
{
  // Iterate along the pairs:
  // - check if everytime at least one of the pair id was listed before
  std::set<uint32_t> previous_ids;
  const Pair_Set used_pairs = Relative_Scale::Get_pairs(relative_scales);
  for (const Pair & pair_it : used_pairs)
  {
    if (previous_ids.empty())
    {
      previous_ids = {pair_it.first, pair_it.second};
    }
    else
      if (previous_ids.count(pair_it.first) == 0 &&
          previous_ids.count(pair_it.second) == 0)
        return false;
  }
  return !relative_scales.empty();
}

/// Solve the relative scale to a common coordinate system
/// Since ratios of depth are used, the found solution is an approximation
/// Solving equation (4) from “Global Structure-from-Motion by Similarity Averaging"
/// Zhaopeng Cui and Ping Tan. (ICCV 2015).”
///
inline
bool
Solve_stellar_translation_scales_averaging
(
  const size_t node_id, // Central node
  const std::vector<Relative_Scale> & vec_relative_scales,
  const Hash_Map<Pair, geometry::Pose3> & relative_poses,
  Hash_Map<IndexT, geometry::Pose3> & triplet_pose,
  const Stellar_Translation_Averaging_Solver_Type e_used_solver =
    Stellar_Translation_Averaging_Solver_Type::SCALING_SOLVER_L2_FULL
)
{
  const Pair_Set used_pairs = Relative_Scale::Get_pairs(vec_relative_scales);
  std::cout << "Stellar reconstruction with center node: " << node_id << "\n"
    << "#relative scales: " << vec_relative_scales.size() << "\n"
    << "#pairs : " << used_pairs.size() << std::endl;

  // Assert that the relative scales pose ids defined a unique connected component
  if (!Relative_scales_are_one_cc(vec_relative_scales))
  {
    std::cerr << "The stellar edges are not giving a single connected component." << std::endl;
    return false;
  }

  // Compute a contiguous indexes mapping
  std::map<Pair, unsigned int> pair_to_index;
  {
    unsigned int i = 0;
    for (const Pair & pair_it: used_pairs)
    {
      pair_to_index[pair_it] = i++;
    }
  }

  // Solve the scale to put all the pair in a common global coordinate system
  Vec x_scales;
  switch (e_used_solver)
  {
    case Stellar_Translation_Averaging_Solver_Type::SCALING_SOLVER_L2_FULL:
    {
      // Setup the linear system: lhx X = rhs
      std::vector< Eigen::Triplet<double> > vec_triplets;
      Vec rhs(vec_relative_scales.size());
      unsigned int i = 0;
      for (const Relative_Scale & relative_scales : vec_relative_scales)
      {
        // Write:
        // S_ij / S_il = factor_ijl
        // as the following linear constraint:
        // S_ij - S_il = log(factor_ijl)
        vec_triplets.emplace_back(i, pair_to_index.at(relative_scales.pairs[0]), 1); // row, col, value
        vec_triplets.emplace_back(i, pair_to_index.at(relative_scales.pairs[1]), -1);
        rhs[i] = std::log(relative_scales.ratio);
        ++i;
      }

      sMat lhs;
      lhs.resize(vec_relative_scales.size(), pair_to_index.size());
      lhs.setFromTriplets(vec_triplets.cbegin(), vec_triplets.cend());

      // Use QR
      Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver(lhs);
      if (solver.info() != Eigen::Success)
      {
        std::cerr << "Sparse matrix cannot be factorized" << std::endl;
        return false;
      }
      const Vec x = solver.solve(rhs);
      if (solver.info() != Eigen::Success)
      {
        std::cerr << "Sparse system cannot be solved" << std::endl;
        return false;
      }

      // Make the log of the distances negative, such that all distances are <= 1
      const double maxLogDistance = x.maxCoeff();
      x_scales = (x.array() - maxLogDistance).array().exp();
    }
    break;

    case Stellar_Translation_Averaging_Solver_Type::SCALING_SOLVER_L2:
    {
      // Setup the linear system: lhx X = rhs
      std::vector< Eigen::Triplet<double> > vec_triplets;
      Vec rhs(vec_relative_scales.size());
      unsigned int i = 0;
      for (const Relative_Scale & relative_scales : vec_relative_scales)
      {
        // Write:
        // S_ij / S_il = factor_ijl
        // as the following linear constraint:
        // S_ij - S_il = log(factor_ijl)
        if (pair_to_index[relative_scales.pairs[0]]!=0)
          vec_triplets.emplace_back(i, pair_to_index.at(relative_scales.pairs[0])-1, 1); // row, col, value
        if (pair_to_index[relative_scales.pairs[1]]!=0)
          vec_triplets.emplace_back(i, pair_to_index.at(relative_scales.pairs[1])-1, -1);
        rhs[i] = std::log(relative_scales.ratio);
        ++i;
      }

      sMat lhs;
      lhs.resize(vec_relative_scales.size(), pair_to_index.size()-1);
      lhs.setFromTriplets(vec_triplets.cbegin(), vec_triplets.cend());

      // Use QR
      Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver(lhs);
      if (solver.info() != Eigen::Success)
      {
        std::cerr << "Sparse matrix cannot be factorized" << std::endl;
        return false;
      }
      const Vec x = solver.solve(rhs);
      if (solver.info() != Eigen::Success)
      {
        std::cerr << "Sparse system cannot be solved" << std::endl;
        return false;
      }

      // Make the log of the distances negative, such that all distances are <= 1
      Vec xCpy(x.size()+1);
      xCpy << 0.0, x;
      const double maxLogDistance = xCpy.maxCoeff();
      x_scales = (xCpy.array() - maxLogDistance).array().exp();
    }
      break;

    case Stellar_Translation_Averaging_Solver_Type::SCALING_SOLVER_L1:
    {
      // Setup the linear system: lhx X = rhs
      std::vector< Eigen::Triplet<double> > vec_triplets;
      Vec rhs(vec_relative_scales.size());
      unsigned int i = 0;
      for (const Relative_Scale & relative_scales : vec_relative_scales)
      {
        // Write:
        // S_ij / S_il = factor_ijl
        // as the following linear constraint:
        // S_ij - S_il = log(factor_ijl)
        if (pair_to_index[relative_scales.pairs[0]]!=0)
          vec_triplets.emplace_back(i, pair_to_index.at(relative_scales.pairs[0])-1, 1); // row, col, value
        if (pair_to_index[relative_scales.pairs[1]]!=0)
          vec_triplets.emplace_back(i, pair_to_index.at(relative_scales.pairs[1])-1, -1);
        rhs[i] = std::log(relative_scales.ratio);
        ++i;
      }

      sMat lhs;
      lhs.resize(vec_relative_scales.size(), pair_to_index.size()-1);
      lhs.setFromTriplets(vec_triplets.cbegin(), vec_triplets.cend());

      // Try to solve under L1 norm
      openMVG::L1Solver<sMat>::Options opt;
      openMVG::L1Solver<sMat> l1solver(opt, lhs);

      Vec l1solution_vec(pair_to_index.size());
      l1solution_vec.setZero();
      if (l1solver.Solve(rhs, &l1solution_vec))
      {
        // Make the log of the distances negative, such that all distances are <= 1
        Vec x(l1solution_vec.size()+1);
        x << 0.0, l1solution_vec;
        const double maxLogDistance = x.maxCoeff();
        x_scales = (x.array() - maxLogDistance).array().exp();
      }
      else
      {
        std::cerr << "Sparse system cannot be solved" << std::endl;
        return false;
      }
    }
    break;
    default:
      std::cerr << "Unsupported solver" << std::endl;
      return false;
  }

  if (x_scales.size() == 0)
  {
    return false;
  }

  // Upgrade the found stellar reconstruction
  triplet_pose.clear();
  triplet_pose[node_id] = geometry::Pose3(); // Identity

  std::set<unsigned int> stellar_ids;
  for (unsigned int i = 0 ; i < x_scales.size(); ++i)
  {
    // Retrieve the original pair (inverse mapping)
    Pair selected_pair;
    for (const Pair & pair_it: used_pairs)
    {
      if (pair_to_index[pair_it] == i)
      {
        selected_pair = pair_it;
        break;
      }
    }
    std::cout << "Pair: " << selected_pair.first << "," << selected_pair.second << "; scaling: " << x_scales[i] << std::endl;
    stellar_ids.insert(selected_pair.first);
    stellar_ids.insert(selected_pair.second);

    if (relative_poses.count(selected_pair) == 0)
    {
      return false;
    }

    if (selected_pair.first == node_id)
    {
      geometry::Pose3 relative_pose = relative_poses.at(selected_pair);
      relative_pose.center() /= x_scales[i];
      triplet_pose[selected_pair.second] = relative_pose;
    }
    else
    {
      geometry::Pose3 relative_pose = relative_poses.at(selected_pair).inverse();
      relative_pose.center() /= x_scales[i];
      triplet_pose[selected_pair.first] = relative_pose;
    }
  }
  return true;
}

} // namespace sfm
} // namespace openMVG


#endif // OPENMVG_SFM_STELLAR_STELLAR_DEFINITIONS_HPP
