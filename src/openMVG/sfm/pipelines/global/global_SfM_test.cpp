
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//-----------------
// Test summary:
//-----------------
// - Create features points and matching from the synthetic dataset
// - Init a SfM_Data scene View and Intrinsic from a synthetic dataset
// - Perform Global SfM on the data
// - Assert that:
//   - mean residual error is below the gaussian noise added to observation
//   - the desired number of tracks are found,
//   - the desired number of poses are found.
//-----------------

#include "openMVG/sfm/pipelines/pipelines_test.hpp"
#include "openMVG/sfm/sfm.hpp"
using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::geometry;
using namespace openMVG::sfm;

#include "testing/testing.h"

#include <cmath>
#include <cstdio>
#include <iostream>

TEST(GLOBAL_SFM, RotationAveragingL2_TranslationAveragingL1) {

  const int nviews = 6;
  const int npoints = 64;
  const nViewDatasetConfigurator config;
  const NViewDataSet d = NRealisticCamerasRing(nviews, npoints, config);

  // Translate the input dataset to a SfM_Data scene
  const SfM_Data sfm_data = getInputScene(d, config, PINHOLE_CAMERA);

  // Remove poses and structure
  SfM_Data sfm_data_2 = sfm_data;
  sfm_data_2.poses.clear();
  sfm_data_2.structure.clear();

  GlobalSfMReconstructionEngine_RelativeMotions sfmEngine(
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

  // Configure motion averaging methods
  sfmEngine.SetRotationAveragingMethod(ROTATION_AVERAGING_L2);
  sfmEngine.SetTranslationAveragingMethod(TRANSLATION_AVERAGING_L1);

  EXPECT_TRUE (sfmEngine.Process());

  const double dResidual = RMSE(sfmEngine.Get_SfM_Data());
  std::cout << "RMSE residual: " << dResidual << std::endl;
  EXPECT_TRUE( dResidual < 0.5);
  EXPECT_EQ( nviews, sfmEngine.Get_SfM_Data().GetPoses().size());
  EXPECT_EQ( npoints, sfmEngine.Get_SfM_Data().GetLandmarks().size());
  EXPECT_TRUE( IsTracksOneCC(sfmEngine.Get_SfM_Data()));
}

TEST(GLOBAL_SFM, RotationAveragingL1_TranslationAveragingL1) {

  const int nviews = 6;
  const int npoints = 64;
  const nViewDatasetConfigurator config;
  const NViewDataSet d = NRealisticCamerasRing(nviews, npoints, config);

  // Translate the input dataset to a SfM_Data scene
  const SfM_Data sfm_data = getInputScene(d, config, PINHOLE_CAMERA);

  // Remove poses and structure
  SfM_Data sfm_data_2 = sfm_data;
  sfm_data_2.poses.clear();
  sfm_data_2.structure.clear();

  GlobalSfMReconstructionEngine_RelativeMotions sfmEngine(
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

  // Configure motion averaging methods
  sfmEngine.SetRotationAveragingMethod(ROTATION_AVERAGING_L1);
  sfmEngine.SetTranslationAveragingMethod(TRANSLATION_AVERAGING_L1);

  EXPECT_TRUE (sfmEngine.Process());

  const double dResidual = RMSE(sfmEngine.Get_SfM_Data());
  std::cout << "RMSE residual: " << dResidual << std::endl;
  EXPECT_TRUE( dResidual < 0.5);
  EXPECT_EQ( nviews, sfmEngine.Get_SfM_Data().GetPoses().size());
  EXPECT_EQ( npoints, sfmEngine.Get_SfM_Data().GetLandmarks().size());
  EXPECT_TRUE( IsTracksOneCC(sfmEngine.Get_SfM_Data()));
}


TEST(GLOBAL_SFM, RotationAveragingL2_TranslationAveragingL2_Chordal) {

  const int nviews = 6;
  const int npoints = 64;
  const nViewDatasetConfigurator config;
  const NViewDataSet d = NRealisticCamerasRing(nviews, npoints, config);

  // Translate the input dataset to a SfM_Data scene
  const SfM_Data sfm_data = getInputScene(d, config, PINHOLE_CAMERA);

  // Remove poses and structure
  SfM_Data sfm_data_2 = sfm_data;
  sfm_data_2.poses.clear();
  sfm_data_2.structure.clear();

  GlobalSfMReconstructionEngine_RelativeMotions sfmEngine(
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

  // Configure motion averaging method
  sfmEngine.SetRotationAveragingMethod(ROTATION_AVERAGING_L2);
  sfmEngine.SetTranslationAveragingMethod(TRANSLATION_AVERAGING_L2_DISTANCE_CHORDAL);

  EXPECT_TRUE (sfmEngine.Process());

  const double dResidual = RMSE(sfmEngine.Get_SfM_Data());
  std::cout << "RMSE residual: " << dResidual << std::endl;
  EXPECT_TRUE( dResidual < 0.5);
  EXPECT_EQ( nviews, sfmEngine.Get_SfM_Data().GetPoses().size());
  EXPECT_EQ( npoints, sfmEngine.Get_SfM_Data().GetLandmarks().size());
  EXPECT_TRUE( IsTracksOneCC(sfmEngine.Get_SfM_Data()));
}

TEST(GLOBAL_SFM, RotationAveragingL2_TranslationAveragingSoftL1) {

  const int nviews = 6;
  const int npoints = 64;
  const nViewDatasetConfigurator config;
  const NViewDataSet d = NRealisticCamerasRing(nviews, npoints, config);

  // Translate the input dataset to a SfM_Data scene
  const SfM_Data sfm_data = getInputScene(d, config, PINHOLE_CAMERA);

  // Remove poses and structure
  SfM_Data sfm_data_2 = sfm_data;
  sfm_data_2.poses.clear();
  sfm_data_2.structure.clear();

  GlobalSfMReconstructionEngine_RelativeMotions sfmEngine(
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

  // Configure motion averaging methods
  sfmEngine.SetRotationAveragingMethod(ROTATION_AVERAGING_L2);
  sfmEngine.SetTranslationAveragingMethod(TRANSLATION_AVERAGING_SOFTL1);

  EXPECT_TRUE (sfmEngine.Process());

  const double dResidual = RMSE(sfmEngine.Get_SfM_Data());
  std::cout << "RMSE residual: " << dResidual << std::endl;
  EXPECT_TRUE( dResidual < 0.5);
  EXPECT_EQ( nviews, sfmEngine.Get_SfM_Data().GetPoses().size());
  EXPECT_EQ( npoints, sfmEngine.Get_SfM_Data().GetLandmarks().size());
  EXPECT_TRUE( IsTracksOneCC(sfmEngine.Get_SfM_Data()));
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
