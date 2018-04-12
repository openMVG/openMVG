// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2018 Yan Qingsong, Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/cameras/Camera_Common.hpp"
#include "openMVG/multiview/test_data_sets.hpp"
#include "openMVG/sfm/sfm.hpp"
#include "openMVG/sfm/pipelines/pipelines_test.hpp"
#include "openMVG/geometry/Similarity3.hpp"
#include "testing/testing.h"
#include <iostream>

using namespace openMVG;
using namespace openMVG::sfm;
using namespace openMVG::cameras;
using namespace openMVG::geometry;
using namespace std;

TEST(SimilarityTest, TestScaleRotTrans) {
  
  // Generate Scene
  const int nviews = 3;
  const int npoints = 10;
  const nViewDatasetConfigurator config;
  const NViewDataSet d = NRealisticCamerasRing(nviews, npoints, config);

  // Translate the input dataset to a SfM_Data scene
  SfM_Data sfm_data = getInputScene(d, config, PINHOLE_CAMERA);
  SfM_Data sfm_data_origin = sfm_data;

  double scale_test;
  Vec3 translation_test, rotation_test;

  // Set random scale,rotaiton and translation
  scale_test = abs(Eigen::MatrixXd::Random(1,1)(0,0))+1;
  rotation_test.setRandom();
  translation_test.setRandom();

  std::cout << "\n"
    << "Scale " << scale_test << "\n"
    << "Rot \n" << rotation_test << "\n"
    << "Tran " << translation_test.transpose();

  Similarity3 sim(Pose3(LookAt(rotation_test),translation_test),scale_test);
  Similarity3 sim_inv = sim.inverse();

  ApplySimilarity(sim,sfm_data);
  ApplySimilarity(sim_inv,sfm_data);
   
  // Check Point
  auto sfm_data_structure_iter = sfm_data.structure.begin();
  auto sfm_data_origin_structure_iter = sfm_data_origin.structure.begin();
  for(;sfm_data_structure_iter!=sfm_data.structure.end() && sfm_data_origin_structure_iter!=sfm_data_origin.structure.end();
       sfm_data_structure_iter++,sfm_data_origin_structure_iter++)
  {
      EXPECT_NEAR(0.0,(sfm_data_structure_iter->second.X-sfm_data_origin_structure_iter->second.X).norm(),1e-8);
  }
  
  // Check Pose
  auto sfm_data_pose_iter = sfm_data.poses.begin();
  auto sfm_data_origin_pose_iter = sfm_data_origin.poses.begin();
  for(;sfm_data_pose_iter!=sfm_data.poses.end() && sfm_data_origin_pose_iter!=sfm_data_origin.poses.end();
       sfm_data_pose_iter++,sfm_data_origin_pose_iter++)
  {
      // Trans
      EXPECT_NEAR(0.0,(sfm_data_pose_iter->second.center()-sfm_data_origin_pose_iter->second.center()).norm(),1e-8);
      // Rot
      EXPECT_NEAR(0.0,(Eigen::AngleAxisd(sfm_data_pose_iter->second.rotation()).matrix()-Eigen::AngleAxisd(sfm_data_origin_pose_iter->second.rotation()).matrix()).norm(),1e-8)
  }
}


/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
