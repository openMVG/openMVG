// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//-----------------
// Test summary:
//-----------------
// - Create features points and matching from the synthetic dataset
// - Init a SfM_Data scene View and Intrinsic from a synthetic dataset
// - Compute some stellar reconstruction on the data
// - Assert that:
//   - The computed poses are correct
//-----------------

#include "openMVG/multiview/essential.hpp"
#include "openMVG/sfm/pipelines/pipelines_test.hpp"
#include "openMVG/sfm/pipelines/relative_pose_engine.hpp"
#include "openMVG/sfm/pipelines/stellar/stellar_solver.hpp"
#include "openMVG/sfm/sfm.hpp"

#include "testing/testing.h"

#include <cmath>
#include <cstdio>
#include <iostream>

using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::geometry;
using namespace openMVG::sfm;

// Test a scene where all the camera intrinsics are known
TEST(STELLAR_SOLVER, Known_Intrinsics) {

  const int nviews = 6;
  const int npoints = 32;
  const nViewDatasetConfigurator config;
  const NViewDataSet d = NRealisticCamerasRing(nviews, npoints, config);

  // Translate the input dataset to a SfM_Data scene
  const SfM_Data sfm_data = getInputScene(d, config, PINHOLE_CAMERA);

  // Remove poses and structure
  SfM_Data sfm_data_2 = sfm_data;
  sfm_data_2.poses.clear();
  sfm_data_2.structure.clear();

  // Configure the features_provider & the matches_provider from the synthetic dataset
  std::shared_ptr<Features_Provider> feats_provider =
    std::make_shared<Synthetic_Features_Provider>();
  // Add a tiny noise in 2D observations to make data more realistic
  std::normal_distribution<double> distribution(0.0, 1e-2);
  dynamic_cast<Synthetic_Features_Provider*>(feats_provider.get())->load(d,distribution);

  std::shared_ptr<Matches_Provider> matches_provider =
    std::make_shared<Synthetic_Matches_Provider>();
  dynamic_cast<Synthetic_Matches_Provider*>(matches_provider.get())->load(d);

  const std::vector<StellarPod> stellar_pods = {
    {{0,1}, {0,2}, {1,2}}, // Test with a triplet
    {{0,1}, {0,2}}         // Test with a 2-uplet
  };

  for (const auto pod : stellar_pods)
  {
    // Compute how many pose ids are defined by the stellar pod
    const std::set<IndexT> set_pose = [&]{
      std::set<IndexT> pose_set;
      for (const Pair & it : pod)
      {
        pose_set.insert(it.first);
        pose_set.insert(it.second);
      }
      return pose_set;
    }();

    const Relative_Pose_Engine::Relative_Pair_Poses relative_poses = [&]
    {
      Relative_Pose_Engine relative_pose_engine;
      if (!relative_pose_engine.Process(
          pod,
          sfm_data_2,
          matches_provider.get(),
          feats_provider.get()))
        return Relative_Pose_Engine::Relative_Pair_Poses();
      else
        return relative_pose_engine.Get_Relative_Poses();
    }();

    EXPECT_EQ(pod.size(), relative_poses.size());

    Stellar_Solver stellar_pod_solver(
      pod,
      relative_poses,
      sfm_data_2,
      matches_provider.get(),
      feats_provider.get());

    EXPECT_TRUE(stellar_pod_solver.Solve(sfm_data_2.poses));

    EXPECT_EQ(set_pose.size(), sfm_data_2.poses.size());

    const double kEpsilon = 1e-3;
    for (const auto & pair : pod)
    {
      const int I = pair.first;
      const int J = pair.second;
      //-- Compute Ground Truth motion
      Mat3 R_gt;
      Vec3 t_gt;
      RelativeCameraMotion(d._R[I], d._t[I], d._R[J], d._t[J], &R_gt, &t_gt);

      // Compute the SfM motion
      Mat3 R_computed;
      Vec3 t_computed;
      RelativeCameraMotion(
        sfm_data_2.poses[I].rotation(),
        sfm_data_2.poses[I].translation(),
        sfm_data_2.poses[J].rotation(),
        sfm_data_2.poses[J].translation(),
        &R_computed, &t_computed);

      // Compare the motion
      std::cout << FrobeniusDistance(R_gt, R_computed) << std::endl;
      std::cout << "translation: " << (t_gt.normalized() - t_computed.normalized()).norm() << std::endl;
      EXPECT_TRUE(FrobeniusDistance(R_gt, R_computed) < kEpsilon);
      EXPECT_TRUE((t_gt.normalized() - t_computed.normalized()).norm() < kEpsilon )
    }
  }
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
