// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2018 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/pipelines/relative_pose_engine.hpp"

#include "openMVG/multiview/essential.hpp"
#include "openMVG/multiview/triangulation.hpp"
#include "openMVG/sfm/pipelines/sfm_robust_model_estimation.hpp"
#include "openMVG/sfm/pipelines/sfm_features_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_matches_provider.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_BA_ceres.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_data_triangulation.hpp"
#include "openMVG/system/timer.hpp"

#include "ceres/ceres.h"

namespace openMVG {
namespace sfm {

using namespace cameras;
using namespace geometry;
using namespace matching;

// Try to compute all the possible relative pose.
bool Relative_Pose_Engine::Process(
  const SfM_Data & sfm_data_,
  const Matches_Provider * matches_provider_,
  const Features_Provider * features_provider_
)
{
  Pair_Set relative_pose_pairs;
  for (const auto & iterMatches : matches_provider_->pairWise_matches_)
  {
    const Pair pair = iterMatches.first;
    const View * v1 = sfm_data_.GetViews().at(pair.first).get();
    const View * v2 = sfm_data_.GetViews().at(pair.second).get();
    if (v1->id_pose != v2->id_pose)
      relative_pose_pairs.insert({v1->id_pose, v2->id_pose});
  }
  return
    Process(
      relative_pose_pairs,
      sfm_data_,
      matches_provider_,
      features_provider_);
}

// Try to compute the relative pose for the provided pose pairs
bool Relative_Pose_Engine::Relative_Pose_Engine::Process(
  const Pair_Set & relative_pose_pairs,
  const SfM_Data & sfm_data_,
  const Matches_Provider * matches_provider_,
  const Features_Provider * features_provider_
)
{
  //
  // List the pairwise matches related to each pose edge ids.
  //
  using PoseWiseMatches = Hash_Map<Pair, Pair_Set>;
  PoseWiseMatches posewise_matches;
  for (const auto & iterMatches : matches_provider_->pairWise_matches_)
  {
    const Pair pair = iterMatches.first;
    const View * v1 = sfm_data_.GetViews().at(pair.first).get();
    const View * v2 = sfm_data_.GetViews().at(pair.second).get();
    if (v1->id_pose != v2->id_pose
        && relative_pose_pairs.count({v1->id_pose, v2->id_pose}) == 1)
      posewise_matches[{v1->id_pose, v2->id_pose}].insert(pair);
  }

  system::Timer t;

  std::unique_ptr<C_Progress> progress_status
    (new C_Progress_display(posewise_matches.size(),
      std::cout, "\n- Relative pose computation -\n" ));

  #ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for schedule(dynamic)
  #endif
  // Compute the relative pose from pairwise point matches:
  for (int i = 0; i < static_cast<int>(posewise_matches.size()); ++i)
  {
    ++(*progress_status);
    {
      PoseWiseMatches::const_iterator iter (posewise_matches.begin());
      std::advance(iter, i);
      const auto & relative_pose_iterator(*iter);
      const Pair relative_pose_pair = relative_pose_iterator.first;
      const Pair_Set & match_pairs = relative_pose_iterator.second;

      // If a pair has the same ID, discard it
      if (relative_pose_pair.first == relative_pose_pair.second)
      {
        continue;
      }

      // Select common bearing vectors
      if (match_pairs.size() > 1)
      {
        std::cerr << "Compute relative pose between more than two view is not supported" << std::endl;
        continue;
      }

      const Pair current_pair(*std::begin(match_pairs));

      const IndexT
        I = current_pair.first,
        J = current_pair.second;

      const View
        * view_I = sfm_data_.views.at(I).get(),
        * view_J = sfm_data_.views.at(J).get();

      // Check that valid cameras exist for the view pair
      if (sfm_data_.GetIntrinsics().count(view_I->id_intrinsic) == 0 ||
          sfm_data_.GetIntrinsics().count(view_J->id_intrinsic) == 0)
        continue;

      const IntrinsicBase
        * cam_I = sfm_data_.GetIntrinsics().at(view_I->id_intrinsic).get(),
        * cam_J = sfm_data_.GetIntrinsics().at(view_J->id_intrinsic).get();

      // Compute for each feature the un-distorted camera coordinates
      const matching::IndMatches & matches = matches_provider_->pairWise_matches_.at(current_pair);
      size_t number_matches = matches.size();
      Mat2X x1(2, number_matches), x2(2, number_matches);
      number_matches = 0;
      for (const auto & match : matches)
      {
        x1.col(number_matches) = cam_I->get_ud_pixel(
          features_provider_->feats_per_view.at(I)[match.i_].coords().cast<double>());
        x2.col(number_matches++) = cam_J->get_ud_pixel(
          features_provider_->feats_per_view.at(J)[match.j_].coords().cast<double>());
      }

      RelativePose_Info relativePose_info;
      relativePose_info.initial_residual_tolerance = Square(2.5);
      if (!robustRelativePose(cam_I, cam_J,
                              x1, x2, relativePose_info,
                              {cam_I->w(), cam_I->h()},
                              {cam_J->w(), cam_J->h()},
                              256))
      {
        continue;
      }
      const bool bRefine_using_BA = true;
      if (bRefine_using_BA)
      {
        // Refine the defined scene
        SfM_Data tiny_scene;
        tiny_scene.views.insert(*sfm_data_.GetViews().find(view_I->id_view));
        tiny_scene.views.insert(*sfm_data_.GetViews().find(view_J->id_view));
        tiny_scene.intrinsics.insert(*sfm_data_.GetIntrinsics().find(view_I->id_intrinsic));
        tiny_scene.intrinsics.insert(*sfm_data_.GetIntrinsics().find(view_J->id_intrinsic));

        // Init poses
        const Pose3 & pose_I = tiny_scene.poses[view_I->id_pose] = {Mat3::Identity(), Vec3::Zero()};
        const Pose3 & pose_J = tiny_scene.poses[view_J->id_pose] = relativePose_info.relativePose;

        // Init structure
        Landmarks & landmarks = tiny_scene.structure;
        for (Mat::Index k = 0; k < x1.cols(); ++k)
        {
          Vec3 X;
          TriangulateDLT(pose_I.asMatrix(), (*cam_I)(x1.col(k)),
                         pose_J.asMatrix(), (*cam_J)(x2.col(k)), &X);
          Observations obs;
          const Vec2 obs_I = features_provider_->feats_per_view.at(I)[matches[k].i_].coords().cast<double>();
          const Vec2 obs_J = features_provider_->feats_per_view.at(J)[matches[k].j_].coords().cast<double>();
          obs[view_I->id_view] = {obs_I, matches[k].i_};
          obs[view_J->id_view] = {obs_J, matches[k].j_};
          landmarks[k].obs = obs;
          landmarks[k].X = X;
        }
        // - refine only Structure and Rotations & translations (keep intrinsic constant)
        Bundle_Adjustment_Ceres::BA_Ceres_options options(false, false);
        options.linear_solver_type_ = ceres::DENSE_SCHUR;
        Bundle_Adjustment_Ceres bundle_adjustment_obj(options);
        const Optimize_Options ba_refine_options
          (Intrinsic_Parameter_Type::NONE, // -> Keep intrinsic constant
          Extrinsic_Parameter_Type::ADJUST_ALL, // adjust camera motion
          Structure_Parameter_Type::ADJUST_ALL);// adjust scene structure
        if (bundle_adjustment_obj.Adjust(tiny_scene, ba_refine_options))
        {
          // --> to debug: save relative pair geometry on disk
          // std::ostringstream os;
          // os << relative_pose_pair.first << "_" << relative_pose_pair.second << ".ply";
          // Save(tiny_scene, os.str(), ESfM_Data(STRUCTURE | EXTRINSICS));
          //
          const Mat3 R1 = tiny_scene.poses[view_I->id_pose].rotation(),
                     R2 = tiny_scene.poses[view_J->id_pose].rotation();
          const Vec3 t1 = tiny_scene.poses[view_I->id_pose].translation(),
                     t2 = tiny_scene.poses[view_J->id_pose].translation();
          // Compute relative motion and save it
          Mat3 Rrel;
          Vec3 trel;
          RelativeCameraMotion(R1, t1, R2, t2, &Rrel, &trel);
          // Update found relative pose
          relativePose_info.relativePose = Pose3(Rrel, -Rrel.transpose() * trel);
        }
      }
#ifdef OPENMVG_USE_OPENMP
      #pragma omp critical
#endif
      {
        // Add the relative pose to the relative 'rotation' pose graph
        relative_poses_[relative_pose_pair] = relativePose_info.relativePose;
      }
    }
  }
  std::cout << "Relative motion computation took: " << t.elapsedMs() << "(ms)" << std::endl;
  return !relative_poses_.empty();
}

// Relative poses accessor
Relative_Pose_Engine::Relative_Pair_Poses
Relative_Pose_Engine::Get_Relative_Poses() const
{
  return relative_poses_;
}

} // namespace sfm
} // namespace openMVG
