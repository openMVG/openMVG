// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/sfm_data_triangulation.hpp"

#include "openMVG/multiview/triangulation_nview.hpp"
#include "openMVG/robust_estimation/rand_sampling.hpp"
#include "third_party/progress/progress.hpp"

#include <deque>
#include <memory>

namespace openMVG {
namespace sfm {

using namespace openMVG::geometry;
using namespace openMVG::cameras;

SfM_Data_Structure_Computation_Basis::SfM_Data_Structure_Computation_Basis
(
  bool bConsoleVerbose
):
  bConsole_verbose_(bConsoleVerbose)
{
}

SfM_Data_Structure_Computation_Blind::SfM_Data_Structure_Computation_Blind
(
  bool bConsoleVerbose
)
  :SfM_Data_Structure_Computation_Basis(bConsoleVerbose)
{
}

void SfM_Data_Structure_Computation_Blind::triangulate
(
  SfM_Data & sfm_data
)
const
{
  std::deque<IndexT> rejectedId;
  std::unique_ptr<C_Progress_display> my_progress_bar;
  if (bConsole_verbose_)
    my_progress_bar.reset( new C_Progress_display(
    sfm_data.structure.size(),
    std::cout,
    "Blind triangulation progress:\n" ));
#ifdef OPENMVG_USE_OPENMP
  #pragma omp parallel
#endif
  for(Landmarks::iterator iterTracks = sfm_data.structure.begin();
    iterTracks != sfm_data.structure.end();
    ++iterTracks)
  {
#ifdef OPENMVG_USE_OPENMP
  #pragma omp single nowait
#endif
    {
      if (bConsole_verbose_)
      {
#ifdef OPENMVG_USE_OPENMP
  #pragma omp critical
#endif
        ++(*my_progress_bar);
      }
      // Triangulate each landmark
      Triangulation trianObj;
      const Observations & obs = iterTracks->second.obs;
      for(Observations::const_iterator itObs = obs.begin();
        itObs != obs.end(); ++itObs)
      {
        const View * view = sfm_data.views.at(itObs->first).get();
        if (sfm_data.IsPoseAndIntrinsicDefined(view))
        {
          const IntrinsicBase * cam = sfm_data.GetIntrinsics().at(view->id_intrinsic).get();
          const Pose3 pose = sfm_data.GetPoseOrDie(view);
          trianObj.add(
            cam->get_projective_equivalent(pose),
            cam->get_ud_pixel(itObs->second.x));
        }
      }
      if (trianObj.size() < 2)
      {
#ifdef OPENMVG_USE_OPENMP
        #pragma omp critical
#endif
        {
          rejectedId.push_front(iterTracks->first);
        }
      }
      else
      {
        // Compute the 3D point
        const Vec3 X = trianObj.compute();
        if (trianObj.minDepth() > 0) // Keep the point only if it have a positive depth
        {
          iterTracks->second.X = X;
        }
        else
        {
#ifdef OPENMVG_USE_OPENMP
          #pragma omp critical
#endif
          {
            rejectedId.push_front(iterTracks->first);
          }
        }
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
  bool bConsoleVerbose
):
  SfM_Data_Structure_Computation_Basis(bConsoleVerbose),
  max_reprojection_error_(max_reprojection_error)
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
    my_progress_bar.reset( new C_Progress_display(
    sfm_data.structure.size(),
    std::cout,
    "Robust triangulation progress:\n" ));
#ifdef OPENMVG_USE_OPENMP
  #pragma omp parallel
#endif
  for(Landmarks::iterator iterTracks = sfm_data.structure.begin();
    iterTracks != sfm_data.structure.end();
    ++iterTracks)
  {
#ifdef OPENMVG_USE_OPENMP
  #pragma omp single nowait
#endif
    {
      if (bConsole_verbose_)
      {
#ifdef OPENMVG_USE_OPENMP
  #pragma omp critical
#endif
        ++(*my_progress_bar);
      }
      Landmark landmark;
      if (robust_triangulation(sfm_data, iterTracks->second.obs, landmark))
      {
        iterTracks->second = landmark;
      }
      else
      {
        // Track must be deleted
#ifdef OPENMVG_USE_OPENMP
  #pragma omp critical
#endif
        {
          rejectedId.push_front(iterTracks->first);
        }
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
/// A point must be seen in at least 3 views
/// Return true for a successful triangulation
bool SfM_Data_Structure_Computation_Robust::robust_triangulation
(
  const SfM_Data & sfm_data,
  const Observations & obs,
  Landmark & landmark, // X & valid observations
  const IndexT min_required_inliers,
  const IndexT min_sample_index
)
const
{
  if (obs.size() < 3)
  {
    return false;
  }

  const double dSquared_pixel_threshold = Square(max_reprojection_error_);

  const IndexT nbIter = obs.size(); // TODO: automatic computation of the number of iterations?

  // - Ransac variables
  Vec3 best_model;
  std::deque<IndexT> best_inlier_set;
  double best_error = std::numeric_limits<double>::max();

  // - Ransac loop
  for (IndexT i = 0; i < nbIter; ++i)
  {
    std::vector<size_t> vec_samples;
    robust::UniformSample(min_sample_index, obs.size(), &vec_samples);
    const std::set<IndexT> samples(vec_samples.begin(), vec_samples.end());

    // Hypothesis generation.
    const Vec3 current_model = track_sample_triangulation(sfm_data, obs, samples);

    // Test validity of the hypothesis
    // - chierality (for the samples)
    // - residual error

    // Chierality (Check the point is in front of the sampled cameras)
    bool bChierality = true;
    bool bReprojection_error = true;
    for (std::set<IndexT>::const_iterator it = samples.begin();
      it != samples.end() && bChierality && bReprojection_error; ++it)
    {
      Observations::const_iterator itObs = obs.begin();
      std::advance(itObs, *it);
      const View * view = sfm_data.views.at(itObs->first).get();
      const IntrinsicBase * cam = sfm_data.GetIntrinsics().at(view->id_intrinsic).get();
      const Pose3 pose = sfm_data.GetPoseOrDie(view);
      const double z = pose.depth(current_model);
      bChierality &= z > 0;
      const Vec2 residual = cam->residual(pose, current_model, itObs->second.x);
      bReprojection_error &= residual.squaredNorm() < dSquared_pixel_threshold;
    }

    if (!bChierality || !bReprojection_error)
      continue;

    std::deque<IndexT> inlier_set;
    double current_error = 0.0;
    // inlier/outlier classification according pixel residual errors.
    for (Observations::const_iterator itObs = obs.begin();
        itObs != obs.end(); ++itObs)
    {
      const View * view = sfm_data.views.at(itObs->first).get();
      const IntrinsicBase * intrinsic = sfm_data.GetIntrinsics().at(view->id_intrinsic).get();
      const Pose3 pose = sfm_data.GetPoseOrDie(view);
      const Vec2 residual = intrinsic->residual(pose, current_model, itObs->second.x);
      const double residual_d = residual.squaredNorm();
      if (residual_d < dSquared_pixel_threshold)
      {
        inlier_set.push_front(itObs->first);
        current_error += residual_d;
      }
      else
      {
        current_error += dSquared_pixel_threshold;
      }
    }
    // Does the hypothesis is the best one we have seen and have sufficient inliers.
    if (current_error < best_error && inlier_set.size() >= min_required_inliers)
    {
      best_model = current_model;
      best_inlier_set = inlier_set;
      best_error = current_error;
    }
  }
  if (!best_inlier_set.empty())
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


/// Triangulate a given track from a selection of observations
Vec3 SfM_Data_Structure_Computation_Robust::track_sample_triangulation(
  const SfM_Data & sfm_data,
  const Observations & obs,
  const std::set<IndexT> & samples) const
{
  Triangulation trianObj;
  for (auto& it : samples)
  {
    const IndexT & idx = it;
    Observations::const_iterator itObs = obs.begin();
    std::advance(itObs, idx);
    const View * view = sfm_data.views.at(itObs->first).get();
    const IntrinsicBase * cam = sfm_data.GetIntrinsics().at(view->id_intrinsic).get();
    const Pose3 pose = sfm_data.GetPoseOrDie(view);
    trianObj.add(
      cam->get_projective_equivalent(pose),
      cam->get_ud_pixel(itObs->second.x));
  }
  return trianObj.compute();
}

} // namespace sfm
} // namespace openMVG
