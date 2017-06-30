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
// - Perform Sequential SfM on the data
// - Assert that:
//   - mean residual error is below the gaussian noise added to observation
//   - the desired number of tracks are found,
//   - the desired number of poses are found.
//-----------------

#include "openMVG/sfm/pipelines/pipelines_test.hpp"
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
TEST(SEQUENTIAL_SFM, Known_Intrinsics) {

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

  SequentialSfMReconstructionEngine sfmEngine(
    sfm_data_2,
    "./",
    stlplus::create_filespec("./", "Reconstruction_Report.html"));

  // Configure the features_provider & the matches_provider from the synthetic dataset
  std::shared_ptr<Features_Provider> feats_provider =
    std::make_shared<Synthetic_Features_Provider>();
  // Add a tiny noise in 2D observations to make data more realistic
  std::normal_distribution<double> distribution(0.0,0.5);
  dynamic_cast<Synthetic_Features_Provider*>(feats_provider.get())->load(d,distribution);

  std::shared_ptr<Matches_Provider> matches_provider =
    std::make_shared<Synthetic_Matches_Provider>();
  dynamic_cast<Synthetic_Matches_Provider*>(matches_provider.get())->load(d);

  // Configure data provider (Features and Matches)
  sfmEngine.SetFeaturesProvider(feats_provider.get());
  sfmEngine.SetMatchesProvider(matches_provider.get());

  // Configure reconstruction parameters (intrinsic parameters are held constant)
  sfmEngine.Set_Intrinsics_Refinement_Type(cameras::Intrinsic_Parameter_Type::NONE);

  // Will use view ids (0,1) as the initial pair
  Views::const_iterator iter_view_0 = sfm_data_2.GetViews().begin();
  Views::const_iterator iter_view_1 = sfm_data_2.GetViews().begin();
  std::advance(iter_view_1, 1);
  sfmEngine.setInitialPair(
    Pair(iter_view_0->second->id_view, iter_view_1->second->id_view));

  EXPECT_TRUE (sfmEngine.Process());

  const double dResidual = RMSE(sfmEngine.Get_SfM_Data());
  std::cout << "RMSE residual: " << dResidual << std::endl;
  EXPECT_TRUE( dResidual < 0.5);
  EXPECT_TRUE( sfmEngine.Get_SfM_Data().GetPoses().size() == nviews);
  EXPECT_TRUE( sfmEngine.Get_SfM_Data().GetLandmarks().size() == npoints);
  EXPECT_TRUE( IsTracksOneCC(sfmEngine.Get_SfM_Data()));
}

// Test a scene where only the two first camera have known intrinsics
TEST(SEQUENTIAL_SFM, Partially_Known_Intrinsics) {

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
  // Only the first three views will have valid intrinsics
  // Remaining one will have undefined intrinsics
  for (Views::iterator iterV = sfm_data_2.views.begin();
    iterV != sfm_data_2.views.end(); ++iterV)
  {
    if (std::distance(sfm_data_2.views.begin(),iterV) > 2)
    {
      iterV->second->id_intrinsic = UndefinedIndexT;
    }
  }

  SequentialSfMReconstructionEngine sfmEngine(
    sfm_data_2,
    "./",
    stlplus::create_filespec("./", "Reconstruction_Report.html"));

  // Configure the features_provider & the matches_provider from the synthetic dataset
  std::shared_ptr<Features_Provider> feats_provider =
    std::make_shared<Synthetic_Features_Provider>();
  // Add a tiny noise in 2D observations to make data more realistic
  std::normal_distribution<double> distribution(0.0,0.5);
  dynamic_cast<Synthetic_Features_Provider*>(feats_provider.get())->load(d,distribution);

  std::shared_ptr<Matches_Provider> matches_provider =
    std::make_shared<Synthetic_Matches_Provider>();
  dynamic_cast<Synthetic_Matches_Provider*>(matches_provider.get())->load(d);

  // Configure data provider (Features and Matches)
  sfmEngine.SetFeaturesProvider(feats_provider.get());
  sfmEngine.SetMatchesProvider(matches_provider.get());

  // Configure reconstruction parameters (intrinsic parameters are held constant)
  sfmEngine.Set_Intrinsics_Refinement_Type(cameras::Intrinsic_Parameter_Type::NONE);

  // Will use view ids (0,2) as the initial pair
  Views::const_iterator iter_view_0 = sfm_data_2.GetViews().begin();
  Views::const_iterator iter_view_2 = sfm_data_2.GetViews().begin();
  std::advance(iter_view_2, 2);
  sfmEngine.setInitialPair(Pair(iter_view_0->second->id_view,
    iter_view_2->second->id_view));

  EXPECT_TRUE (sfmEngine.Process());

  const double dResidual = RMSE(sfmEngine.Get_SfM_Data());
  std::cout << "RMSE residual: " << dResidual << std::endl;
  EXPECT_TRUE( dResidual < 0.55);
  EXPECT_TRUE( sfmEngine.Get_SfM_Data().GetPoses().size() == nviews);
  EXPECT_TRUE( sfmEngine.Get_SfM_Data().GetLandmarks().size() == npoints);
  EXPECT_TRUE( IsTracksOneCC(sfmEngine.Get_SfM_Data()));
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
