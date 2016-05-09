
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/multiview/test_data_sets.hpp"
#include "openMVG/sfm/sfm.hpp"
using namespace openMVG;
using namespace openMVG::sfm;

#include <random>
#include <iostream>

// Create from a synthetic scene (NViewDataSet) some SfM pipelines data provider:
//  - for each view store the observations point as PointFeatures
struct Synthetic_Features_Provider : public Features_Provider
{
  template <typename NoiseGenerator>
  bool load(
    const NViewDataSet & synthetic_data, NoiseGenerator & noise)
  {
    std::default_random_engine generator;
    // For each view
    for (int j = 0; j < synthetic_data._n; ++j)
    {
      // For each new point visibility
      for (int i = 0; i < synthetic_data._x[j].cols(); ++i)
      {
        const Vec2 pt = synthetic_data._x[j].col(i);
        feats_per_view[j].push_back(
          features::PointFeature(pt(0)+noise(generator), pt(1)+noise(generator)));
      }
    }
    return true;
  }
};

// Create from a synthetic scene (NViewDataSet) some SfM pipelines data provider:
//  - for contiguous triplets store the corresponding observations indexes
struct Synthetic_Matches_Provider : public Matches_Provider
{
  virtual bool load(
    const NViewDataSet & synthetic_data)
  {
    // For each view
    for (int j = 0; j < synthetic_data._n; ++j)
    {
      for (int jj = j+1; jj < j+3 ; ++jj)
      {
        for (int idx = 0; idx < synthetic_data._x[j].cols(); ++idx)
        {
          pairWise_matches_[Pair(j,(jj)%synthetic_data._n)].push_back(IndMatch(idx,idx));
        }
      }
    }
    return true;
  }
};

/// Compute the Root Mean Square Error of the residuals
static double RMSE(const SfM_Data & sfm_data)
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
      const geometry::Pose3 pose = sfm_data.GetPoseOrDie(view);
      const std::shared_ptr<cameras::IntrinsicBase> intrinsic = sfm_data.GetIntrinsics().at(view->id_intrinsic);
      const Vec2 residual = intrinsic->residual(pose, iterTracks->second.X, itObs->second.x);
      //std::cout << residual << " ";
      vec.push_back( residual(0) );
      vec.push_back( residual(1) );
    }
  }
  const Eigen::Map<Eigen::RowVectorXd> residuals(&vec[0], vec.size());
  const double RMSE = std::sqrt(residuals.squaredNorm() / vec.size());
  return RMSE;
}

// Translate a synthetic scene into a valid SfM_Data scene
// As only one intrinsic is defined we used shared intrinsic
SfM_Data getInputScene
(
  const NViewDataSet & d,
  const nViewDatasetConfigurator & config,
  cameras::EINTRINSIC eintrinsic)
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
    sfm_data.views[i] = std::make_shared<View>
      ("", id_view, id_intrinsic, id_pose, config._cx *2, config._cy *2);
  }

  // 2. Poses
  for (int i = 0; i < nviews; ++i)
  {
    sfm_data.poses[i] = geometry::Pose3(d._R[i], d._C[i]);;
  }

  // 3. Intrinsic data (shared, so only one camera intrinsic is defined)
  {
    const unsigned int w = config._cx *2;
    const unsigned int h = config._cy *2;
    switch (eintrinsic)
    {
      case cameras::PINHOLE_CAMERA:
        sfm_data.intrinsics[0] = std::make_shared<cameras::Pinhole_Intrinsic>
          (w, h, config._fx, config._cx, config._cy);
      break;
      case cameras::PINHOLE_CAMERA_RADIAL1:
        sfm_data.intrinsics[0] = std::make_shared<cameras::Pinhole_Intrinsic_Radial_K1>
          (w, h, config._fx, config._cx, config._cy, 0.0);
      break;
      case cameras::PINHOLE_CAMERA_RADIAL3:
        sfm_data.intrinsics[0] = std::make_shared<cameras::Pinhole_Intrinsic_Radial_K3>
          (w, h, config._fx, config._cx, config._cy, 0., 0., 0.);
      break;
      default:
        std::cout << "Not yet supported" << std::endl;
    }
  }

  // 4. Landmarks
  for (int i = 0; i < npoints; ++i) {
    // Collect the image of point i in each frame.
    Landmark landmark;
    landmark.X = d._X.col(i);
    for (int j = 0; j < nviews; ++j) {
      const Vec2 pt = d._x[j].col(i);
      landmark.obs[j] = Observation(pt, i);
    }
    sfm_data.structure[i] = landmark;
  }

  return sfm_data;
}

