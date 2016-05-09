
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

#include "openMVG/multiview/test_data_sets.hpp"
#include "openMVG/sfm/sfm.hpp"
using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::geometry;
using namespace openMVG::sfm;

#include "testing/testing.h"
#include "../cameras/Camera_Common.hpp"

#include <cmath>
#include <cstdio>
#include <iostream>

double RMSE(const SfM_Data & sfm_data);

SfM_Data getInputScene(const NViewDataSet & d, const nViewDatasetConfigurator & config, EINTRINSIC eintrinsic);

TEST(BUNDLE_ADJUSTMENT, EffectiveMinimization_Pinhole) {

  const int nviews = 3;
  const int npoints = 6;
  const nViewDatasetConfigurator config;
  const NViewDataSet d = NRealisticCamerasRing(nviews, npoints, config);

  // Translate the input dataset to a SfM_Data scene
  SfM_Data sfm_data = getInputScene(d, config, PINHOLE_CAMERA);

  const double dResidual_before = RMSE(sfm_data);

  // Call the BA interface and let it refine (Structure and Camera parameters [Intrinsics|Motion])
  std::shared_ptr<Bundle_Adjustment> ba_object = std::make_shared<Bundle_Adjustment_Ceres>();
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
  std::shared_ptr<Bundle_Adjustment> ba_object = std::make_shared<Bundle_Adjustment_Ceres>();
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
  std::shared_ptr<Bundle_Adjustment> ba_object = std::make_shared<Bundle_Adjustment_Ceres>();
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
  std::shared_ptr<Bundle_Adjustment> ba_object = std::make_shared<Bundle_Adjustment_Ceres>();
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
  std::shared_ptr<Bundle_Adjustment> ba_object = std::make_shared<Bundle_Adjustment_Ceres>();
  EXPECT_TRUE( ba_object->Adjust(sfm_data,
    Optimize_Options(
      Intrinsic_Parameter_Type::ADJUST_ALL,
      Extrinsic_Parameter_Type::ADJUST_ALL,
      Structure_Parameter_Type::ADJUST_ALL)) );

  const double dResidual_after = RMSE(sfm_data);
  EXPECT_TRUE( dResidual_before > dResidual_after);
}

/// Compute the Root Mean Square Error of the residuals
double RMSE(const SfM_Data & sfm_data)
{
  // Compute residuals for each observation
  std::vector<double> vec;
  for(Landmarks::const_iterator iterTracks = sfm_data.GetLandmarks().begin();
      iterTracks != sfm_data.GetLandmarks().end();
      ++iterTracks)
  {
    const Observations & obs = iterTracks->second.obs;
    for(Observations::const_iterator itObs = obs.begin();
      itObs != obs.end(); ++itObs)
    {
      const View * view = sfm_data.GetViews().find(itObs->first)->second.get();
      const Pose3 pose = sfm_data.GetPoseOrDie(view);
      const std::shared_ptr<IntrinsicBase> intrinsic = sfm_data.GetIntrinsics().find(view->id_intrinsic)->second;
      const Vec2 residual = intrinsic->residual(pose, iterTracks->second.X, itObs->second.x);
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
//    a random noise between [-.5,.5] is added on observed data points
SfM_Data getInputScene(const NViewDataSet & d, const nViewDatasetConfigurator & config, EINTRINSIC eintrinsic)
{
  // Translate the input dataset to a SfM_Data scene
  SfM_Data sfm_data;

  // 1. Views
  // 2. Poses
  // 3. Intrinsic data (shared, so only one camera intrinsic is defined)
  // 4. Landmarks

  const int nviews = d._C.size();
  const int npoints = d._X.cols();

  // 1. Views
  for (int i = 0; i < nviews; ++i)
  {
    const IndexT id_view = i, id_pose = i, id_intrinsic = 0; //(shared intrinsics)
    sfm_data.views[i] = std::make_shared<View>("", id_view, id_intrinsic, id_pose, config._cx *2, config._cy *2);
  }

  // 2. Poses
  for (int i = 0; i < nviews; ++i)
  {
    Pose3 pose(d._R[i], d._C[i]);
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
  for (int i = 0; i < npoints; ++i) {
    // Collect observation of the landmarks X in each frame.
    Landmark landmark;
    landmark.X = d._X.col(i);
    for (int j = 0; j < nviews; ++j) {
      Vec2 pt = d._x[j].col(i);
      // => random noise between [-.5,.5] is added
      pt(0) += rand()/RAND_MAX - .5;
      pt(1) += rand()/RAND_MAX - .5;

      landmark.obs[j] = Observation(pt, i);
    }
    sfm_data.structure[i] = landmark;
  }

  return sfm_data;
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
