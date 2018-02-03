// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/sfm_data_triangulation.hpp"

#include <deque>

#include "openMVG/geometry/pose3.hpp"
#include "openMVG/multiview/triangulation_nview.hpp"
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

/// Triangulate a given track from a selection of observations
bool track_sample_triangulation
(
  const SfM_Data & sfm_data,
  const Observations & obs,
  const std::set<IndexT> & samples,
  Vec3 & X
)
{
  if (samples.size() >= 2 && obs.size() >= 2)
  {
    std::vector<Vec3> bearing;
    std::vector<Mat34> poses;
    bearing.reserve(samples.size());
    poses.reserve(samples.size());
    for (const auto& idx : samples)
    {
      Observations::const_iterator itObs = obs.begin();
      std::advance(itObs, idx);
      const View * view = sfm_data.views.at(itObs->first).get();
      if (!sfm_data.IsPoseAndIntrinsicDefined(view))
        continue;
      const IntrinsicBase * cam = sfm_data.GetIntrinsics().at(view->id_intrinsic).get();
      const Pose3 pose = sfm_data.GetPoseOrDie(view);
      bearing.emplace_back((*cam)(cam->get_ud_pixel(itObs->second.x)));
      poses.emplace_back(pose.asMatrix());
    }
    if (bearing.size() >= 2)
    {
      const Eigen::Map<const Mat3X> bearing_matrix(bearing[0].data(), 3, bearing.size());
      Vec4 Xhomogeneous;
      TriangulateNViewAlgebraic
      (
        bearing_matrix,
        poses, // Ps are projective cameras.
        &Xhomogeneous);
      X = Xhomogeneous.hnormalized();
      return true;
    }
  }
  return false;
}


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
        std::set<IndexT> samples;
        for (size_t i = 0; i < obs.size(); ++i)  { samples.insert(i); }
        Vec3 X;
        if (track_sample_triangulation(sfm_data, obs, samples, X))
        {
          bool bChierality = true;
          for (Observations::const_iterator obs_it = obs.begin();
            obs_it != obs.end() && bChierality; ++obs_it)
          {
            const View * view = sfm_data.views.at(obs_it->first).get();
            if (!sfm_data.IsPoseAndIntrinsicDefined(view))
              continue;
            const IntrinsicBase * cam = sfm_data.intrinsics.at(view->id_intrinsic).get();
            const Pose3 pose = sfm_data.GetPoseOrDie(view);
            bChierality &= CheiralityTest((*cam)(obs_it->second.x), pose, X);
          }

          if (bChierality) // Keep the point only if it has a positive depth
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
  bool bConsoleVerbose
):
  SfM_Data_Structure_Computation_Basis(bConsoleVerbose),
  max_reprojection_error_(max_reprojection_error),
  min_required_inliers_(min_required_inliers),
  min_sample_index_(min_sample_index)
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

/// Robustly try to estimate the best 3D point using a ransac Scheme
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

  // Handle the case where all observations must be used
  if (min_required_inliers_ == min_sample_index_ &&
      obs.size() == min_required_inliers_)
  {
    std::set<IndexT> samples;
    for (int i = 0; i < min_required_inliers_; ++i)  { samples.insert(i); }
    // Generate the 3D point hypothesis by triangulating the observations
    Vec3 X;
    if (track_sample_triangulation(sfm_data, obs, samples, X))
    {
      // Test validity of the hypothesis:
      // - residual error
      // - chierality
      bool bChierality = true;
      bool bReprojection_error = true;
      IndexT validity_test_count = 0;
      for (std::set<IndexT>::const_iterator it = samples.begin();
        it != samples.end() && bChierality && bReprojection_error; ++it)
      {
        Observations::const_iterator itObs = obs.begin();
        std::advance(itObs, *it);
        const View * view = sfm_data.views.at(itObs->first).get();
        if (!sfm_data.IsPoseAndIntrinsicDefined(view))
          continue;
        const IntrinsicBase * cam = sfm_data.GetIntrinsics().at(view->id_intrinsic).get();
        const Pose3 pose = sfm_data.GetPoseOrDie(view);
        bChierality &= CheiralityTest((*cam)(itObs->second.x), pose, X);
        const Vec2 residual = cam->residual(pose(X), itObs->second.x);
        bReprojection_error &= residual.squaredNorm() < dSquared_pixel_threshold;
        validity_test_count += (bChierality && bReprojection_error) ? 1 : 0;
      }

      if (bChierality && bReprojection_error &&
          validity_test_count >= min_required_inliers_)
      {
        landmark.X = X;
        landmark.obs = obs;
        return true;
      }
    }
    return false;
  }

  // We must perform a robust estimation
  // - There is more observations than the minimal number of required sample

  const IndexT nbIter = obs.size(); // TODO: automatic computation of the number of iterations?

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
    std::vector<uint32_t> vec_samples;
    robust::UniformSample(min_sample_index_, obs.size(), random_generator, &vec_samples);
    const std::set<IndexT> samples(vec_samples.begin(), vec_samples.end());

    // Hypothesis generation
    Vec3 X;
    if (!track_sample_triangulation(sfm_data, obs, samples, X))
      continue;

    // Test validity of the hypothesis
    // - chierality (for the samples)
    // - residual error

    bool bChierality = true;
    bool bReprojection_error = true;
    IndexT validity_test_count = 0;
    for (std::set<IndexT>::const_iterator it = samples.begin();
      it != samples.end() && bChierality && bReprojection_error; ++it)
    {
      Observations::const_iterator itObs = obs.begin();
      std::advance(itObs, *it);
      const View * view = sfm_data.views.at(itObs->first).get();
      if (!sfm_data.IsPoseAndIntrinsicDefined(view))
        continue;
      const IntrinsicBase * cam = sfm_data.GetIntrinsics().at(view->id_intrinsic).get();
      const Pose3 pose = sfm_data.GetPoseOrDie(view);
      bChierality &= CheiralityTest((*cam)(itObs->second.x), pose, X);;
      const Vec2 residual = cam->residual(pose(X), itObs->second.x);
      bReprojection_error &= residual.squaredNorm() < dSquared_pixel_threshold;
      validity_test_count += (bChierality && bReprojection_error) ? 1 : 0;
    }

    if (!bChierality || !bReprojection_error ||
        validity_test_count < min_required_inliers_)
      continue;

    std::deque<IndexT> inlier_set;
    double current_error = 0.0;
    // inlier/outlier classification according pixel residual errors.
    for (const auto & obs_it : obs)
    {
      const View * view = sfm_data.views.at(obs_it.first).get();
      if (!sfm_data.IsPoseAndIntrinsicDefined(view))
        continue;
      const IntrinsicBase * intrinsic = sfm_data.GetIntrinsics().at(view->id_intrinsic).get();
      const Pose3 pose = sfm_data.GetPoseOrDie(view);
      const Vec2 residual = intrinsic->residual(pose(X), obs_it.second.x);
      const double residual_d = residual.squaredNorm();
      if (residual_d < dSquared_pixel_threshold)
      {
        inlier_set.push_front(obs_it.first);
        current_error += residual_d;
      }
      else
      {
        current_error += dSquared_pixel_threshold;
      }
    }
    // Does the hypothesis is the best one we have seen and have sufficient inliers.
    if (current_error < best_error && inlier_set.size() >= min_required_inliers_)
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
