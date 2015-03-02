
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


// An example of a minimal, self-contained bundle adjuster using Ceres
// It refines Structure & Camera parameters [Intrinsic|Rotation|Translation].
// => A synthetic scene is used:
//    a random noise between [-.5,.5] is added on observed data points

#include "testing/testing.h"

#include "openMVG/multiview/test_data_sets.hpp"
#include "openMVG/multiview/projection.hpp"

#include "openMVG/sfm/sfm.hpp"

using namespace openMVG;

#include <cmath>
#include <cstdio>
#include <iostream>

double RMSE(const SfM_Data & sfm_data);

TEST(BUNDLE_ADJUSTMENT, EffectiveMinimization_RTf) {

  const int nviews = 3;
  const int npoints = 6;
  const nViewDatasetConfigurator config;
  const NViewDataSet d = NRealisticCamerasRing(nviews, npoints, config);

  // Translate the input dataset to a SfM_Data scene
  SfM_Data sfm_data;

  // 1. Views
  // 2. Poses
  // 3. Intrinsic data (shared, so only one camera intrinsic is defined)
  // 4. Landmarks

  // 1. Views
  for (int i = 0; i < nviews; ++i)
  {
    View v;
    v.id_intrinsic = 0;
    v.id_view = v.id_pose = i;
    sfm_data.views[i] = v;
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
    sfm_data.intrinsics[0] = std::make_shared<Intrinsic>(w, h, config._fx, config._cx, config._cy);
  }

  // 4. Landmarks
  for (int i = 0; i < npoints; ++i) {
    // Collect the image of point i in each frame.
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
  const double dResidual_before = RMSE(sfm_data);

  // Call the BA interface and let it refine (Structure and Camera parameters [Intrinsics|Motion])
  std::shared_ptr<Bundle_Adjustment> ba_object = std::make_shared<Bundle_Adjustment_Ceres>();
  EXPECT_TRUE( ba_object->bAdjust(sfm_data) );

  const double dResidual_after = RMSE(sfm_data);
  EXPECT_TRUE( dResidual_before > dResidual_after);

}

/// Projection of a 3D point into the camera plane (Apply pose, disto and Intrinsics)
static Vec2 Project(
  const geometry::Pose3 & pose,
  const IntrinsicBase * intrinsic,
  const Vec3 & pt3D)
{
  const Vec3 X = pose(pt3D); // apply pose
  if (intrinsic->have_disto()) // apply disto & intrinsics
    return intrinsic->cam2ima( intrinsic->apply(X.head<2>()/X(2)) );
  else // apply intrinsics
    return intrinsic->cam2ima( X.head<2>()/X(2) );
}

/// Compute the Root Mean Square Error of the residuals
double RMSE(const SfM_Data & sfm_data)
{
  // Compute residuals for each observation
  IndexT index = 0;
  std::vector<double> vec;
  for(Landmarks::const_iterator iterTracks = sfm_data.getLandmarks().begin();
      iterTracks != sfm_data.getLandmarks().end();
      ++iterTracks)
  {
    const Observations & obs = iterTracks->second.obs;
    for(Observations::const_iterator itObs = obs.begin();
      itObs != obs.end(); ++itObs, ++index)
    {
      const View & view = sfm_data.getViews().find(itObs->first)->second;
      const Pose3 & pose = sfm_data.getPoses().find(view.id_pose)->second;
      const std::shared_ptr<IntrinsicBase> intrinsic = sfm_data.getIntrinsics().find(view.id_intrinsic)->second;
      const Vec2 proj = Project(pose, intrinsic.get(), iterTracks->second.X);
      const Vec2 residual = proj - itObs->second.x;
      vec.push_back( residual(0) );
      vec.push_back( residual(1) );
    }
  }
  const Eigen::Map<Eigen::RowVectorXd> residuals(&vec[0], vec.size());
  const double RMSE = std::sqrt(residuals.squaredNorm() / vec.size());
  return RMSE;
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
