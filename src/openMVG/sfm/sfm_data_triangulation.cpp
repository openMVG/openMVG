// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/sfm_data_triangulation.hpp"

#include <deque>
#include <functional>

#include "openMVG/geometry/pose3.hpp"
#include "openMVG/multiview/triangulation_nview.hpp"
#include "openMVG/multiview/triangulation.hpp"
#include "openMVG/robust_estimation/rand_sampling.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_landmark.hpp"

#include "third_party/progress/progress_display.hpp"

namespace openMVG {
namespace sfm {

using namespace openMVG::geometry;
using namespace openMVG::cameras;

SfM_Data_Structure_Computation_Basis::SfM_Data_Structure_Computation_Basis
(
  bool bConsoleVerbose
)
  :bConsole_verbose_(bConsoleVerbose)
{
}

SfM_Data_Structure_Computation_Blind::SfM_Data_Structure_Computation_Blind
(
  bool bConsoleVerbose
)
  :SfM_Data_Structure_Computation_Basis(bConsoleVerbose)
{
}

/// Triangulate a given set of observations
bool track_triangulation
(
  const SfM_Data & sfm_data,
  const Observations & obs,
  Vec3 & X,
  const ETriangulationMethod & etri_method = ETriangulationMethod::DEFAULT
)
{
  if (obs.size() >= 2)
  {
    std::vector<Vec3> bearing;
    std::vector<Mat34> poses;
    std::vector<Pose3> poses_;
    bearing.reserve(obs.size());
    poses.reserve(obs.size());
    for (const auto& observation : obs)
    {
      const View * view = sfm_data.views.at(observation.first).get();
      if (!sfm_data.IsPoseAndIntrinsicDefined(view))
        return false;
      const IntrinsicBase * cam = sfm_data.GetIntrinsics().at(view->id_intrinsic).get();
      const Pose3 pose = sfm_data.GetPoseOrDie(view);
      bearing.emplace_back((*cam)(cam->get_ud_pixel(observation.second.x)));
      poses.emplace_back(pose.asMatrix());
      poses_.emplace_back(pose);
    }
    if (bearing.size() > 2)
    {
      const Eigen::Map<const Mat3X> bearing_matrix(bearing[0].data(), 3, bearing.size());
      Vec4 Xhomogeneous;
      if (TriangulateNViewAlgebraic
      (
        bearing_matrix,
        poses,
        &Xhomogeneous))
      {
        X = Xhomogeneous.hnormalized();
        return true;
      }
    }
    else
    {
      return Triangulate2View
      (
        poses_.front().rotation(),
        poses_.front().translation(),
        bearing.front(),
        poses_.back().rotation(),
        poses_.back().translation(),
        bearing.back(),
        X,
        etri_method
      );
    }
  }
  return false;
}

// Test if a predicate is true for each observation
// i.e: predicate could be:
// - cheirality test (depth test): cheirality_predicate
// - cheirality and residual error: ResidualAndCheiralityPredicate::predicate
bool track_check_predicate
(
  const Observations & obs,
  const SfM_Data & sfm_data,
  const Vec3 & X,
  std::function<bool(
    const IntrinsicBase&,
    const Pose3&,
    const Vec2&,
    const Vec3&)> predicate
)
{
  bool visibility = false; // assume that no observation has been looked yet
  for (const auto & obs_it : obs)
  {
    const View * view = sfm_data.views.at(obs_it.first).get();
    if (!sfm_data.IsPoseAndIntrinsicDefined(view))
      continue;
    visibility = true; // at least an observation is evaluated
    const IntrinsicBase * cam = sfm_data.intrinsics.at(view->id_intrinsic).get();
    const Pose3 pose = sfm_data.GetPoseOrDie(view);
    if (!predicate(*cam, pose, obs_it.second.x, X))
      return false;
  }
  return visibility;
}

bool cheirality_predicate
(
  const IntrinsicBase& cam,
  const Pose3& pose,
  const Vec2& x,
  const Vec3& X
)
{
  return CheiralityTest(cam(x), pose, X);
}

struct ResidualAndCheiralityPredicate
{
  const double squared_pixel_threshold_;

  ResidualAndCheiralityPredicate(const double squared_pixel_threshold)
    :squared_pixel_threshold_(squared_pixel_threshold){}

  bool predicate
  (
    const IntrinsicBase& cam,
    const Pose3& pose,
    const Vec2& x,
    const Vec3& X
  )
  {
    const Vec2 residual = cam.residual(pose(X), x);
    return CheiralityTest(cam(x), pose, X) &&
           residual.squaredNorm() < squared_pixel_threshold_;
  }
};

void SfM_Data_Structure_Computation_Blind::triangulate
(
  SfM_Data & sfm_data
)
const
{
  std::deque<IndexT> rejectedId;
  std::unique_ptr<C_Progress> my_progress_bar;
  if (bConsole_verbose_)
    my_progress_bar.reset(
      new C_Progress_display(
        sfm_data.structure.size(),
        std::cout,
        "Blind triangulation progress:\n" ));
#ifdef OPENMVG_USE_OPENMP
  #pragma omp parallel
#endif
  for (auto& tracks_it :sfm_data.structure)
  {
#ifdef OPENMVG_USE_OPENMP
  #pragma omp single nowait
#endif
    {
      if (bConsole_verbose_)
      {
        ++(*my_progress_bar);
      }

      const Observations & obs = tracks_it.second.obs;
      bool bKeep = false;
      {
        // Generate the track 3D hypothesis
        Vec3 X;
        if (track_triangulation(sfm_data, obs, X))
        {
          // Keep the point only if it has a positive depth for all obs
          if (track_check_predicate(obs, sfm_data, X, cheirality_predicate))
          {
            tracks_it.second.X = X;
            bKeep = true;
          }
        }
      }
      if (!bKeep)
      {
#ifdef OPENMVG_USE_OPENMP
        #pragma omp critical
#endif
        rejectedId.push_front(tracks_it.first);
      }
    }
  }
  // Erase the unsuccessful triangulated tracks
  for (auto& it : rejectedId)
  {
    sfm_data.structure.erase(it);
  }
}

SfM_Data_Structure_Computation_Robust::SfM_Data_Structure_Computation_Robust
(
  const double max_reprojection_error,
  const IndexT min_required_inliers,
  const IndexT min_sample_index,
  const ETriangulationMethod etri_method,
  bool bConsoleVerbose
):
  SfM_Data_Structure_Computation_Basis(bConsoleVerbose),
  max_reprojection_error_(max_reprojection_error),
  min_required_inliers_(min_required_inliers),
  min_sample_index_(min_sample_index),
  etri_method_(etri_method)
{
}

void SfM_Data_Structure_Computation_Robust::triangulate
(
  SfM_Data & sfm_data
)
const
{
  robust_triangulation(sfm_data);
}

/// Robust triangulation of track data contained in the structure
/// All observations must have View with valid Intrinsic and Pose data
/// Invalid landmark are removed.
void SfM_Data_Structure_Computation_Robust::robust_triangulation
(
  SfM_Data & sfm_data
)
const
{
  std::deque<IndexT> rejectedId;
  std::unique_ptr<C_Progress_display> my_progress_bar;
  if (bConsole_verbose_)
    my_progress_bar.reset(
      new C_Progress_display(
        sfm_data.structure.size(),
        std::cout,
        "Robust triangulation progress:\n" ));
#ifdef OPENMVG_USE_OPENMP
  #pragma omp parallel
#endif
  for (auto& tracks_it :sfm_data.structure)
  {
#ifdef OPENMVG_USE_OPENMP
  #pragma omp single nowait
#endif
    {
      if (bConsole_verbose_)
      {
        ++(*my_progress_bar);
      }
      Landmark landmark;
      if (robust_triangulation(sfm_data, tracks_it.second.obs, landmark))
      {
        tracks_it.second = landmark;
      }
      else
      {
        // Track must be deleted
#ifdef OPENMVG_USE_OPENMP
        #pragma omp critical
#endif
        rejectedId.push_front(tracks_it.first);
      }
    }
  }
  // Erase the unsuccessful triangulated tracks
  for (auto& it : rejectedId)
  {
    sfm_data.structure.erase(it);
  }
}

Observations ObservationsSampler
(
  const Observations & obs,
  const std::vector<std::uint32_t> & samples
)
{
  Observations sampled_obs;
  for (const auto& idx : samples)
  {
    Observations::const_iterator obs_it = obs.cbegin();
    std::advance(obs_it, idx);
    sampled_obs.insert(*obs_it);
  }
  return sampled_obs;
}

/// Robustly try to estimate the best 3D point using a ransac scheme
/// A point must be seen in at least min_required_inliers views
/// Return true for a successful triangulation
bool SfM_Data_Structure_Computation_Robust::robust_triangulation
(
  const SfM_Data & sfm_data,
  const Observations & obs,
  Landmark & landmark // X & valid observations
)
const
{
  if (obs.size() < min_required_inliers_ || obs.size() < min_sample_index_)
  {
    return false;
  }

  const double dSquared_pixel_threshold = Square(max_reprojection_error_);

  // Predicate to validate a sample (cheirality and residual error)
  ResidualAndCheiralityPredicate predicate(dSquared_pixel_threshold);
  auto predicate_binding = std::bind(&ResidualAndCheiralityPredicate::predicate,
                                     predicate,
                                     std::placeholders::_1,
                                     std::placeholders::_2,
                                     std::placeholders::_3,
                                     std::placeholders::_4);

  // Handle the case where all observations must be used
  if (min_required_inliers_ == min_sample_index_ &&
      obs.size() == min_required_inliers_)
  {
    // Generate the 3D point hypothesis by triangulating all the observations
    Vec3 X;
    if (track_triangulation(sfm_data, obs, X, etri_method_) &&
        track_check_predicate(obs, sfm_data, X, predicate_binding))
    {
      landmark.X = X;
      landmark.obs = obs;
      return true;
    }
    return false;
  }

  // else we perform a robust estimation since
  //  there is more observations than the minimal number of required sample.

  const IndexT nbIter = obs.size() * 2; // TODO: automatic computation of the number of iterations?

  // - Ransac variables
  Vec3 best_model = Vec3::Zero();
  std::deque<IndexT> best_inlier_set;
  double best_error = std::numeric_limits<double>::max();

  //--
  // Random number generation
  std::mt19937 random_generator(std::mt19937::default_seed);

  // - Ransac loop
  for (IndexT i = 0; i < nbIter; ++i)
  {
    std::vector<uint32_t> samples;
    robust::UniformSample(min_sample_index_, obs.size(), random_generator, &samples);

    Vec3 X;
    // Hypothesis generation
    const auto minimal_sample = ObservationsSampler(obs, samples);

    if (!track_triangulation(sfm_data, minimal_sample, X, etri_method_))
      continue;

    // Test validity of the hypothesis
    if (!track_check_predicate(minimal_sample, sfm_data, X, predicate_binding))
      continue;

    std::deque<IndexT> inlier_set;
    double current_error = 0.0;
    // inlier/outlier classification according pixel residual errors.
    for (const auto & obs_it : obs)
    {
      const View * view = sfm_data.views.at(obs_it.first).get();
      if (!sfm_data.IsPoseAndIntrinsicDefined(view))
        continue;
      const IntrinsicBase & cam =  *sfm_data.GetIntrinsics().at(view->id_intrinsic).get();
      const Pose3 pose = sfm_data.GetPoseOrDie(view);
      if (!CheiralityTest(cam(obs_it.second.x), pose, X))
        continue;
      const double residual_sq = cam.residual(pose(X), obs_it.second.x).squaredNorm();
      if (residual_sq < dSquared_pixel_threshold)
      {
        inlier_set.push_front(obs_it.first);
        current_error += residual_sq;
      }
      else
      {
        current_error += dSquared_pixel_threshold;
      }
    }
    // Does the hypothesis:
    // - is the best one we have seen so far.
    // - has sufficient inliers.
    if (current_error < best_error &&
      inlier_set.size() >= min_required_inliers_)
    {
      best_model = X;
      best_inlier_set = inlier_set;
      best_error = current_error;
    }
  }
  if (!best_inlier_set.empty() && best_inlier_set.size() >= min_required_inliers_)
  {
    // Update information (3D landmark position & valid observations)
    landmark.X = best_model;
    for (const IndexT & val : best_inlier_set)
    {
      landmark.obs[val] = obs.at(val);
    }
  }
  return !best_inlier_set.empty();
}

} // namespace sfm
} // namespace openMVG
