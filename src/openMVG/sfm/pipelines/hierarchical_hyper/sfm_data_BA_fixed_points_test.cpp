// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 nomoko AG, Sebastien Chappuis<sebastien@nomoko.camera>, Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/multiview/test_data_sets.hpp"
#include "openMVG/sfm/pipelines/hierarchical_hyper/sfm_data_BA_fixed_points.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/cameras/cameras.hpp"

#include "testing/testing.h"

using namespace openMVG;
using namespace openMVG::sfm;
using namespace openMVG::cameras;
using namespace std;

SfM_Data getInputScene(const NViewDataSet & data, const nViewDatasetConfigurator & config, EINTRINSIC eintrinsic);
double RMSE(const SfM_Data & sfm_data);

TEST(Bundle_Adjustment_Fixed_Points, effectiveMinimization)
{
  // create a scene
  const int nviews = 3;
  const int npoints = 6;
  const nViewDatasetConfigurator config;
  const NViewDataSet d = NRealisticCamerasRing(nviews, npoints, config);

  // Convert the input dataset to a SfM_Data scene
  SfM_Data sfm_data = getInputScene(d, config, PINHOLE_CAMERA_RADIAL3);

  const double dResidual_before = RMSE(sfm_data);

  const std::set<IndexT> fixed_tracks_ids{0,1};

  Bundle_Adjustment_Fixed_Points bundle_adjustment_obj;

  bundle_adjustment_obj.Adjust(sfm_data, fixed_tracks_ids);

  const double dResidual_after = RMSE(sfm_data);
  EXPECT_TRUE( dResidual_before > dResidual_after);
}

TEST(Bundle_Adjustment_Fixed_Points, leavesFixedTracksConstant)
{
  // create a scene
  const int nviews = 3;
  const int npoints = 6;
  const nViewDatasetConfigurator config;
  const NViewDataSet d = NRealisticCamerasRing(nviews, npoints, config);

  // Convert the input dataset to a SfM_Data scene
  SfM_Data sfm_data = getInputScene(d, config, PINHOLE_CAMERA_RADIAL3);

  const std::set<IndexT> fixed_tracks_ids = {0,1};

  const Landmarks landmarks_before = sfm_data.GetLandmarks();

  Bundle_Adjustment_Fixed_Points bundle_adjustment_obj;

  bundle_adjustment_obj.Adjust(sfm_data, fixed_tracks_ids);

  const Landmarks landmarks_after = sfm_data.GetLandmarks();

  // check that fixed landmarks stayed fixed and that the others "moved"
  EXPECT_EQ(landmarks_before.at(0).X, landmarks_after.at(0).X);
  EXPECT_EQ(landmarks_before.at(1).X, landmarks_after.at(1).X);
  EXPECT_FALSE((landmarks_before.at(2).X == landmarks_after.at(2).X));
}

// NOTE : STOLEN FROM src/openMVG/sfm/sfm_data_BA_test.cpp
/// Compute the Root Mean Square Error of the residuals
double RMSE(const SfM_Data & sfm_data)
{
  // Compute residuals for each observation
  std::vector<double> vec;
  for (Landmarks::const_iterator iterTracks = sfm_data.GetLandmarks().begin();
      iterTracks != sfm_data.GetLandmarks().end();
      ++iterTracks)
  {
    const Observations & obs = iterTracks->second.obs;
    for (Observations::const_iterator itObs = obs.begin();
      itObs != obs.end(); ++itObs)
    {
      const View * view = sfm_data.GetViews().find(itObs->first)->second.get();
      const Pose3 pose = sfm_data.GetPoseOrDie(view);
      const std::shared_ptr<IntrinsicBase> intrinsic = sfm_data.GetIntrinsics().find(view->id_intrinsic)->second;
      const Vec2 residual = intrinsic->residual(pose(iterTracks->second.X), itObs->second.x);
      vec.push_back( residual(0) );
      vec.push_back( residual(1) );
    }
  }
  const Eigen::Map<Eigen::RowVectorXd> residuals(&vec[0], vec.size());
  const double RMSE = std::sqrt(residuals.squaredNorm() / vec.size());
  return RMSE;
}

// NOTE : STOLEN FROM src/openMVG/sfm/sfm_data_BA_test.cpp
// Translation a synthetic scene into a valid SfM_Data scene.
// => A synthetic scene is used:
//    a random noise between [-.5,.5] is added on observed data points
SfM_Data getInputScene(const NViewDataSet & data, const nViewDatasetConfigurator & config, EINTRINSIC eintrinsic)
{
  // Translate the input dataset to a SfM_Data scene
  SfM_Data sfm_data;

  // 1. Views
  // 2. Poses
  // 3. Intrinsic data (shared, so only one camera intrinsic is defined)
  // 4. Landmarks

  const int nviews = data._C.size();
  const int npoints = data._X.cols();

  // 1. Views
  for (int i = 0; i < nviews; ++i)
  {
    const IndexT id_view = i, id_pose = i, id_intrinsic = 0; //(shared intrinsics)
    sfm_data.views[i] = std::make_shared<View>("", id_view, id_intrinsic, id_pose, config._cx *2, config._cy *2);
  }

  // 2. Poses
  for (int i = 0; i < nviews; ++i)
  {
    Pose3 pose(data._R[i], data._C[i]);
    sfm_data.poses[i] = pose;
  }

  // 3. Intrinsic data
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
    landmark.X = data._X.col(i);
    for (int j = 0; j < nviews; ++j) {
      Vec2 pt = data._x[j].col(i);
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
