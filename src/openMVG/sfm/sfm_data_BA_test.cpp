// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//-----------------
// Test summary:
//-----------------
// - Create a SfM_Data scene from a synthetic dataset
//   - since random noise have been added on 2d data point (initial residual is not small)
// - Check that residual is small once the generic Bundle Adjustment framework have been called.
// --
// - Perform the test for all the plausible intrinsic camera models
//-----------------

#include "openMVG/cameras/Camera_Common.hpp"
#include "openMVG/multiview/test_data_sets.hpp"
#include "openMVG/sfm/sfm.hpp"

#include "testing/testing.h"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <random>

using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::geometry;
using namespace openMVG::sfm;

double RMSE(const SfM_Data & sfm_data);

SfM_Data getInputScene
(
  const NViewDataSet & d,
  const nViewDatasetConfigurator & config,
  EINTRINSIC eintrinsic,
  const bool b_use_gcp = false,
  const bool b_use_pose_prior = false,
  const bool b_use_noise_on_image_observations = true
);

TEST(BUNDLE_ADJUSTMENT, EffectiveMinimization_Pinhole) {

  const int nviews = 3;
  const int npoints = 6;
  const nViewDatasetConfigurator config;
  const NViewDataSet d = NRealisticCamerasRing(nviews, npoints, config);

  // Translate the input dataset to a SfM_Data scene
  SfM_Data sfm_data = getInputScene(d, config, PINHOLE_CAMERA);

  const double dResidual_before = RMSE(sfm_data);

  // Call the BA interface and let it refine (Structure and Camera parameters [Intrinsics|Motion])
  const bool bVerbose = true;
  const bool bMultithread = false;
  std::shared_ptr<Bundle_Adjustment> ba_object =
    std::make_shared<Bundle_Adjustment_Ceres>(
      Bundle_Adjustment_Ceres::BA_Ceres_options(bVerbose, bMultithread));
  EXPECT_TRUE( ba_object->Adjust(sfm_data,
    Optimize_Options(
      Intrinsic_Parameter_Type::ADJUST_ALL,
      Extrinsic_Parameter_Type::ADJUST_ALL,
      Structure_Parameter_Type::ADJUST_ALL)) );

  const double dResidual_after = RMSE(sfm_data);
  EXPECT_TRUE( dResidual_before > dResidual_after);
}

TEST(BUNDLE_ADJUSTMENT, EffectiveMinimization_Pinhole_Radial_K1) {

  const int nviews = 3;
  const int npoints = 6;
  const nViewDatasetConfigurator config;
  const NViewDataSet d = NRealisticCamerasRing(nviews, npoints, config);

  // Translate the input dataset to a SfM_Data scene
  SfM_Data sfm_data = getInputScene(d, config, PINHOLE_CAMERA_RADIAL1);

  const double dResidual_before = RMSE(sfm_data);

  // Call the BA interface and let it refine (Structure and Camera parameters [Intrinsics|Motion])
  const bool bVerbose = true;
  const bool bMultithread = false;
  std::shared_ptr<Bundle_Adjustment> ba_object =
    std::make_shared<Bundle_Adjustment_Ceres>(
      Bundle_Adjustment_Ceres::BA_Ceres_options(bVerbose, bMultithread));
  EXPECT_TRUE( ba_object->Adjust(sfm_data,
    Optimize_Options(
      Intrinsic_Parameter_Type::ADJUST_ALL,
      Extrinsic_Parameter_Type::ADJUST_ALL,
      Structure_Parameter_Type::ADJUST_ALL)) );

  const double dResidual_after = RMSE(sfm_data);
  EXPECT_TRUE( dResidual_before > dResidual_after);
}

TEST(BUNDLE_ADJUSTMENT, EffectiveMinimization_Pinhole_Radial_K3) {

  const int nviews = 3;
  const int npoints = 6;
  const nViewDatasetConfigurator config;
  const NViewDataSet d = NRealisticCamerasRing(nviews, npoints, config);

  // Translate the input dataset to a SfM_Data scene
  SfM_Data sfm_data = getInputScene(d, config, PINHOLE_CAMERA_RADIAL3);

  const double dResidual_before = RMSE(sfm_data);

  // Call the BA interface and let it refine (Structure and Camera parameters [Intrinsics|Motion])
  const bool bVerbose = true;
  const bool bMultithread = false;
  std::shared_ptr<Bundle_Adjustment> ba_object =
    std::make_shared<Bundle_Adjustment_Ceres>(
      Bundle_Adjustment_Ceres::BA_Ceres_options(bVerbose, bMultithread));
  EXPECT_TRUE( ba_object->Adjust(sfm_data,
    Optimize_Options(
      Intrinsic_Parameter_Type::ADJUST_ALL,
      Extrinsic_Parameter_Type::ADJUST_ALL,
      Structure_Parameter_Type::ADJUST_ALL)) );

  const double dResidual_after = RMSE(sfm_data);
  EXPECT_TRUE( dResidual_before > dResidual_after);
}

TEST(BUNDLE_ADJUSTMENT, EffectiveMinimization_Pinhole_Intrinsic_Brown_T2) {

  const int nviews = 3;
  const int npoints = 6;
  const nViewDatasetConfigurator config;
  const NViewDataSet d = NRealisticCamerasRing(nviews, npoints, config);

  // Translate the input dataset to a SfM_Data scene
  SfM_Data sfm_data = getInputScene(d, config, PINHOLE_CAMERA_BROWN);

  const double dResidual_before = RMSE(sfm_data);

  // Call the BA interface and let it refine (Structure and Camera parameters [Intrinsics|Motion])
  const bool bVerbose = true;
  const bool bMultithread = false;
  std::shared_ptr<Bundle_Adjustment> ba_object =
    std::make_shared<Bundle_Adjustment_Ceres>(
      Bundle_Adjustment_Ceres::BA_Ceres_options(bVerbose, bMultithread));
  EXPECT_TRUE( ba_object->Adjust(sfm_data,
    Optimize_Options(
      Intrinsic_Parameter_Type::ADJUST_ALL,
      Extrinsic_Parameter_Type::ADJUST_ALL,
      Structure_Parameter_Type::ADJUST_ALL)) );

  const double dResidual_after = RMSE(sfm_data);
  EXPECT_TRUE( dResidual_before > dResidual_after);
}

TEST(BUNDLE_ADJUSTMENT, EffectiveMinimization_Pinhole_Intrinsic_Fisheye) {

  const int nviews = 3;
  const int npoints = 6;
  const nViewDatasetConfigurator config;
  const NViewDataSet d = NRealisticCamerasRing(nviews, npoints, config);

  // Translate the input dataset to a SfM_Data scene
  SfM_Data sfm_data = getInputScene(d, config, PINHOLE_CAMERA_FISHEYE);

  const double dResidual_before = RMSE(sfm_data);

  // Call the BA interface and let it refine (Structure and Camera parameters [Intrinsics|Motion])
  const bool bVerbose = true;
  const bool bMultithread = false;
  std::shared_ptr<Bundle_Adjustment> ba_object =
    std::make_shared<Bundle_Adjustment_Ceres>(
      Bundle_Adjustment_Ceres::BA_Ceres_options(bVerbose, bMultithread));
  EXPECT_TRUE( ba_object->Adjust(sfm_data,
    Optimize_Options(
      Intrinsic_Parameter_Type::ADJUST_ALL,
      Extrinsic_Parameter_Type::ADJUST_ALL,
      Structure_Parameter_Type::ADJUST_ALL)) );

  const double dResidual_after = RMSE(sfm_data);
  EXPECT_TRUE( dResidual_before > dResidual_after);
}

//-- Test with GCP - Camera position once BA done must be the same as the GT
TEST(BUNDLE_ADJUSTMENT, EffectiveMinimization_Pinhole_GCP) {

  const int nviews = 3;
  const int npoints = 6;
  const nViewDatasetConfigurator config;
  const NViewDataSet d = NRealisticCamerasRing(nviews, npoints, config);

  // Translate the input dataset to a SfM_Data scene
  SfM_Data sfm_data = getInputScene(d, config, PINHOLE_CAMERA, true);

  const double dResidual_before = RMSE(sfm_data);

  // Call the BA interface and let it refine (Structure and Camera parameters [Intrinsics|Motion])
  const bool bVerbose = true;
  const bool bMultithread = false;
  std::shared_ptr<Bundle_Adjustment> ba_object =
    std::make_shared<Bundle_Adjustment_Ceres>(
      Bundle_Adjustment_Ceres::BA_Ceres_options(bVerbose, bMultithread));
  EXPECT_TRUE( ba_object->Adjust(sfm_data,
    Optimize_Options(
      Intrinsic_Parameter_Type::NONE,
      Extrinsic_Parameter_Type::ADJUST_ALL,
      Structure_Parameter_Type::ADJUST_ALL,
      //-> Use GCP to fit the SfM scene to the GT coordinates system (Scale, Rotation, Translation)
      Control_Point_Parameter(20.0, true))) );

  const double dResidual_after = RMSE(sfm_data);
  EXPECT_TRUE( dResidual_before > dResidual_after);

  //-- Check camera pose are to the right place (since GCP was used, the camera coordinates must be the same)
  for (const auto & view_it : sfm_data.GetViews())
  {
    const View * view = view_it.second.get();
    const Pose3 pose = sfm_data.GetPoseOrDie(view);
    const double position_residual = (d._C[view->id_pose] - pose.center()).norm();
    EXPECT_NEAR(0.0, position_residual, 1e-4);
  }
}

//-- Test with PosePriors - Camera position once BA done must be the same as the GT
TEST(BUNDLE_ADJUSTMENT, EffectiveMinimization_Pinhole_PosePriors) {

  const int nviews = 12;
  const int npoints = 6;
  const nViewDatasetConfigurator config;
  const NViewDataSet d = NRealisticCamerasRing(nviews, npoints, config);

  // Translate the input dataset to a SfM_Data scene
  const bool b_use_GCP = false;
  const bool b_use_POSE_PRIOR = true;
  const bool b_use_noise_on_image_observations = false;
  SfM_Data sfm_data = getInputScene(d, config, PINHOLE_CAMERA,
    b_use_GCP, b_use_POSE_PRIOR, b_use_noise_on_image_observations);

  const double dResidual_before = RMSE(sfm_data);

  // First run a BA without pose prior
  // - check that RMSE is tiny and that the scene is not at the right place (the pose prior positions)

  {
    const bool b_use_POSE_PRIOR_in_BA = false;

    const bool bVerbose = true;
    const bool bMultithread = false;
    std::shared_ptr<Bundle_Adjustment> ba_object =
      std::make_shared<Bundle_Adjustment_Ceres>(
        Bundle_Adjustment_Ceres::BA_Ceres_options(bVerbose, bMultithread));
    EXPECT_TRUE( ba_object->Adjust(sfm_data,
      Optimize_Options(
        Intrinsic_Parameter_Type::ADJUST_ALL,
        Extrinsic_Parameter_Type::ADJUST_ALL,
        Structure_Parameter_Type::ADJUST_ALL,
        Control_Point_Parameter(0.0, false),
        b_use_POSE_PRIOR_in_BA)) );

    const double dResidual_after = RMSE(sfm_data);
    EXPECT_TRUE( dResidual_before > dResidual_after);

    // Compute distance between SfM poses center and GPS Priors
    Mat residuals(3, sfm_data.GetViews().size());
    for (const auto & view_it : sfm_data.GetViews())
    {
      const ViewPriors * view = dynamic_cast<ViewPriors*>(view_it.second.get());
      residuals.col(view->id_pose) = (sfm_data.GetPoseOrDie(view).center() - view->pose_center_).transpose();
    }
    // Check that the scene is not at the position of the POSE PRIORS center.
    EXPECT_FALSE( (residuals.colwise().norm().sum() < 1e-8) );
  }

  // Then activate BA with pose prior & check that the scene is at the right place
  {
    const bool b_use_POSE_PRIOR_in_BA = true;

    const bool bVerbose = true;
    const bool bMultithread = false;
    std::shared_ptr<Bundle_Adjustment> ba_object =
      std::make_shared<Bundle_Adjustment_Ceres>(
        Bundle_Adjustment_Ceres::BA_Ceres_options(bVerbose, bMultithread));
    EXPECT_TRUE( ba_object->Adjust(sfm_data,
      Optimize_Options(
        Intrinsic_Parameter_Type::ADJUST_ALL,
        Extrinsic_Parameter_Type::ADJUST_ALL,
        Structure_Parameter_Type::ADJUST_ALL,
        Control_Point_Parameter(0.0, false),
        b_use_POSE_PRIOR_in_BA)) );

    const double dResidual_after = RMSE(sfm_data);
    EXPECT_TRUE( dResidual_before > dResidual_after);

    // Compute distance between SfM poses center and GPS Priors
    for (const auto & view_it : sfm_data.GetViews())
    {
      const ViewPriors * view = dynamic_cast<ViewPriors*>(view_it.second.get());
      const Pose3 pose = sfm_data.GetPoseOrDie(view);
      const double position_residual = (d._C[view->id_pose] - pose.center()).norm();
      EXPECT_NEAR(0.0, position_residual, 1e-8);
    }
  }
}


/// Compute the Root Mean Square Error of the residuals
double RMSE(const SfM_Data & sfm_data)
{
  // Compute residuals for each observation
  std::vector<double> vec;
  for (const auto& landmark_it : sfm_data.GetLandmarks())
  {
    const Observations & obs = landmark_it.second.obs;
    for (const auto& obs_it : obs)
    {
      const View * view = sfm_data.GetViews().find(obs_it.first)->second.get();
      const Pose3 pose = sfm_data.GetPoseOrDie(view);
      const std::shared_ptr<IntrinsicBase> intrinsic = sfm_data.GetIntrinsics().find(view->id_intrinsic)->second;
      const Vec2 residual = intrinsic->residual(pose(landmark_it.second.X), obs_it.second.x);
      vec.push_back( residual(0) );
      vec.push_back( residual(1) );
    }
  }
  const Eigen::Map<Eigen::RowVectorXd> residuals(&vec[0], vec.size());
  const double RMSE = std::sqrt(residuals.squaredNorm() / vec.size());
  return RMSE;
}

// Translation a synthetic scene into a valid SfM_Data scene.
// => A synthetic scene is used:
//    some random noise is added on observed structure data points
//    a tiny rotation to ground truth is added to the true rotation (in order to test BA effectiveness)
SfM_Data getInputScene
(
  const NViewDataSet & d,
  const nViewDatasetConfigurator & config,
  EINTRINSIC eintrinsic,
  const bool b_use_gcp,
  const bool b_use_pose_prior,
  const bool b_use_noise_on_image_observations
)
{
  // Translate the input dataset to a SfM_Data scene
  SfM_Data sfm_data;

  // 1. Views
  // 2. Poses
  // 3. Intrinsic data (shared, so only one camera intrinsic is defined)
  // 4. Landmarks
  // 5. GCP (optional)

  const int nviews = d._C.size();
  const int npoints = d._X.cols();

  // 1. Views
  for (int i = 0; i < nviews; ++i)
  {
    const IndexT id_view = i, id_pose = i, id_intrinsic = 0; //(shared intrinsics)

    if (!b_use_pose_prior)
    {
      sfm_data.views[i] = std::make_shared<View>("", id_view, id_intrinsic, id_pose, config._cx *2, config._cy *2);
    }
    else // b_use_pose_prior == true
    {
      sfm_data.views[i] = std::make_shared<ViewPriors>("", id_view, id_intrinsic, id_pose, config._cx *2, config._cy *2);
      ViewPriors * view = dynamic_cast<ViewPriors*>(sfm_data.views[i].get());
      view->b_use_pose_center_ = true;
      view->pose_center_ = d._C[i];
    }
  }

  // Add a rotation to the GT (in order to make BA do some work)
  const Mat3 rot = RotationAroundX(D2R(6));

  // 2. Poses
  for (int i = 0; i < nviews; ++i)
  {
    const Pose3 pose(rot * d._R[i], d._C[i]);
    sfm_data.poses[i] = pose;
  }

  // 3. Intrinsic data (shared, so only one camera intrinsic is defined)
  {
    const unsigned int w = config._cx *2;
    const unsigned int h = config._cy *2;
    switch (eintrinsic)
    {
      case PINHOLE_CAMERA:
        sfm_data.intrinsics[0] = std::make_shared<Pinhole_Intrinsic>
          (w, h, config._fx, config._cx, config._cy);
      break;
      case PINHOLE_CAMERA_RADIAL1:
        sfm_data.intrinsics[0] = std::make_shared<Pinhole_Intrinsic_Radial_K1>
          (w, h, config._fx, config._cx, config._cy, 0.0);
      break;
      case PINHOLE_CAMERA_RADIAL3:
        sfm_data.intrinsics[0] = std::make_shared<Pinhole_Intrinsic_Radial_K3>
          (w, h, config._fx, config._cx, config._cy, 0., 0., 0.);
      break;
      case PINHOLE_CAMERA_BROWN:
        sfm_data.intrinsics[0] = std::make_shared<Pinhole_Intrinsic_Brown_T2>
          (w, h, config._fx, config._cx, config._cy, 0., 0., 0., 0., 0.);
      break;
      case PINHOLE_CAMERA_FISHEYE:
      sfm_data.intrinsics[0] = std::make_shared<Pinhole_Intrinsic_Fisheye>
          (w, h, config._fx, config._cx, config._cy, 0., 0., 0., 0.);
      break;
      default:
        std::cout << "Not yet supported" << std::endl;
    }
  }

  // 4. Landmarks
  // Collect image observation of the landmarks X in each frame.
  // => add some random noise to each (x,y) observation
  std::default_random_engine random_generator;
  std::normal_distribution<double> distribution(0, 0.1);
  for (int i = 0; i < npoints; ++i)
  {
    // Create a landmark for each 3D points
    Landmark landmark;
    landmark.X = d._X.col(i);
    for (int j = 0; j < nviews; ++j)
    {
      Vec2 pt = d._x[j].col(i);
      if (b_use_noise_on_image_observations)
      {
        // Add some noise to image observations
        pt(0) += distribution(random_generator);
        pt(1) += distribution(random_generator);
      }

      landmark.obs[j] = Observation(pt, i);
    }
    sfm_data.structure[i] = landmark;
  }

  // 5. GCP
  if (b_use_gcp)
  {
    if (npoints >= 4) // Use 4 GCP for this test
    {
      for (int i = 0; i < 4; ++i) // Select the 4 first point as GCP
      {
        // Collect observations of the landmarks X in each frame.
        Landmark landmark;
        landmark.X = d._X.col(i);
        for (int j = 0; j < nviews; ++j)
        {
          landmark.obs[j] = Observation(d._x[j].col(i), i);
        }
        sfm_data.control_points[i] = landmark;
      }
    }
    else
    {
      std::cerr << "Insufficient point count" << std::endl;
    }
  }
  return sfm_data;
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
