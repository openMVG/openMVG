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

#include "openMVG/sfm/base/pipelines_test.hpp"
#include "openMVG/sfm/sfm.hpp"

#include "openMVG/multiview/trifocal/trifocal_model.hpp"
#include "openMVG/multiview/trifocal/solver_trifocal_metrics.hpp"
#include "testing/testing.h"

#include <cmath>
#include <cstdio>
#include <iostream>

using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::geometry;
using namespace openMVG::sfm;

static void
invert_intrinsics(
    const double K[/*3 or 2 ignoring last line*/][3],
    const double px_coords[2],
    double normalized_coords[2])
{
  const double *px = px_coords;
  double *nrm = normalized_coords;
  nrm[1] = (px[1] - K[1][2]) /K[1][1];
  nrm[0] = (px[0] - K[0][1]*nrm[1] - K[0][2])/K[0][0];
}

//static void
//invert_intrinsics_tgt(
//    const double K[/*3 or 2 ignoring last line*/][3],
//    const double px_tgt_coords[2],
//    double normalized_tgt_coords[2])
//{
//  const double *tp = px_tgt_coords;
//  double *t = normalized_tgt_coords;
//  t[1] = tp[1]/K[1][1];
//  t[0] = (tp[0] - K[0][1]*t[1])/K[0][0];
//  double n = hypot(t[0], t[1]);
//  t[0] /= n; t[1] /= n;
//}
#if 0
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
  sfmEngine.setInitialPair({sfm_data_2.GetViews().at(0)->id_view,
                            sfm_data_2.GetViews().at(1)->id_view});

  EXPECT_TRUE (sfmEngine.Process());

  const double dResidual = RMSE(sfmEngine.Get_SfM_Data());
  std::cout << "RMSE residual: " << dResidual << std::endl;
  EXPECT_TRUE( dResidual < 0.5);
  EXPECT_TRUE( sfmEngine.Get_SfM_Data().GetPoses().size() == nviews);
  EXPECT_TRUE( sfmEngine.Get_SfM_Data().GetLandmarks().size() == npoints);
  EXPECT_TRUE( IsTracksOneCC(sfmEngine.Get_SfM_Data()));
}
#endif

#if 0
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
  for (auto & view : sfm_data_2.views)
  {
    if (view.second->id_view > 2)
    {
      view.second->id_intrinsic = UndefinedIndexT;
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
  std::normal_distribution<double> distribution(0.0, 0.5);
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
  sfmEngine.setInitialPair({sfm_data_2.GetViews().at(0)->id_view,
                            sfm_data_2.GetViews().at(2)->id_view});

  EXPECT_TRUE (sfmEngine.Process());

  const double dResidual = RMSE(sfmEngine.Get_SfM_Data());
  std::cout << "RMSE residual: " << dResidual << std::endl;
  EXPECT_TRUE( dResidual < 0.55);
  EXPECT_TRUE( sfmEngine.Get_SfM_Data().GetPoses().size() == nviews);
  EXPECT_TRUE( sfmEngine.Get_SfM_Data().GetLandmarks().size() == npoints);
  EXPECT_TRUE( IsTracksOneCC(sfmEngine.Get_SfM_Data()));
}
#endif

#if 1
//- Oriented tests ------------------------------------------------------------

//-----------------
// Test summary:
//-----------------
// - Create oriented features points and matching from the synthcurves synthetic dataset
//    - 5 features x 5 true point correspondences x 4 views
// - Init a SfM_Data scene View and Intrinsic from a synthetic dataset
// - Perform Sequential SfM on the data
// - Assert that:
//   - mean residual error is below the gaussian noise added to observation
//   - the desired number of tracks are found,
//   - the desired number of poses are found.
//-----------------

// Create from a synthetic scene (NViewDataSet) some SfM pipelines data provider:
//  - for each view store the observations point as PointFeatures
struct Synthetic_Oriented_Features_Provider : public Features_Provider
{
  bool load( const NViewOrientedDataSet & synthetic_data) {
    // For each view
    for (size_t v = 0; v < synthetic_data._n; ++v) {
      // For each new point visibility
      for (Mat2X::Index i = 0; i < synthetic_data._x[v].cols(); ++i) {
        const Vec2 pt = synthetic_data._x[v].col(i);
        const Vec2 tgt = synthetic_data._tgt2d[v].col(i);
        sio_feats_per_view[v].emplace_back(pt(0), pt(1), 1.0, atan2(tgt(1),tgt(0)));
      }
    }
    return true;
  }
};


bool
check_camera_triplet(const NViewOrientedDataSet &d, const int ci[]) 
{
  const int nviews = d.nviews();
  const int npoints = d.npts();

  openMVG::trifocal::trifocal_model_t gt_cam;
  std::vector<Vec4> datum0(npoints);
  std::vector<Vec4> datum1(npoints);
  std::vector<Vec4> datum2(npoints);
  Mat3 K;
  // Fill data
  for(unsigned r = 0; r < 3; ++r)
    for(unsigned c = 0; c < 3; ++c) {
      gt_cam[0](r,c) = d._cameras_gt_raw[ci[0]][r][c];
      gt_cam[1](r,c) = d._cameras_gt_raw[ci[1]][r][c];
      gt_cam[2](r,c) = d._cameras_gt_raw[ci[2]][r][c];
      if (r != 3)
        K(r,c) = d._K_raw[r][c]; 
    }
  K(2,1) = K(2,0) = 0;
  K(2,2) = 1;
  for(unsigned r = 0; r < 3; ++r) {
    gt_cam[0](r,3) = -(d._cameras_gt_raw[ci[0]][r][0]*d._cameras_gt_raw[ci[0]][3][0]+d._cameras_gt_raw[ci[0]][r][1]*d._cameras_gt_raw[ci[0]][3][1]+d._cameras_gt_raw[ci[0]][r][2]*d._cameras_gt_raw[ci[0]][3][2]);
    gt_cam[1](r,3) = -(d._cameras_gt_raw[ci[1]][r][0]*d._cameras_gt_raw[ci[1]][3][0]+d._cameras_gt_raw[ci[1]][r][1]*d._cameras_gt_raw[ci[1]][3][1]+d._cameras_gt_raw[ci[1]][r][2]*d._cameras_gt_raw[ci[1]][3][2]);
    gt_cam[2](r,3) = -(d._cameras_gt_raw[ci[2]][r][0]*d._cameras_gt_raw[ci[2]][3][0]+d._cameras_gt_raw[ci[2]][r][1]*d._cameras_gt_raw[ci[2]][3][1]+d._cameras_gt_raw[ci[2]][r][2]*d._cameras_gt_raw[ci[2]][3][2]);
  }
  for (unsigned v = 0; v < 3; v++)
    OPENMVG_LOG_INFO << "\n "<< gt_cam[v];

  for (unsigned p = 0; p < npoints; ++p) {
    for (unsigned i = 0; i < 2; ++i) { // Separate pt from tgt
      datum0[p](i)   = d._x[ci[0]].col(p)(i);
      datum0[p](i+2) = d._tgt2d[ci[0]].col(p)(i);
      datum1[p](i)   = d._x[ci[1]].col(p)(i);
      datum1[p](i+2) = d._tgt2d[ci[1]].col(p)(i);
      datum2[p](i)   = d._x[ci[2]].col(p)(i);
      datum2[p](i+2) = d._tgt2d[ci[2]].col(p)(i);
    }
    invert_intrinsics(d._K_raw, datum0[p].data(), datum0[p].data()); 
    invert_intrinsics_tgt(d._K_raw, datum0[p].data()+2, datum0[p].data()+2); 
    invert_intrinsics(d._K_raw, datum1[p].data(), datum1[p].data()); 
    invert_intrinsics_tgt(d._K_raw, datum1[p].data()+2, datum1[p].data()+2); 
    invert_intrinsics(d._K_raw, datum2[p].data(), datum2[p].data()); 
    invert_intrinsics_tgt(d._K_raw, datum2[p].data()+2, datum2[p].data()+2);
    datum0[p].tail(2) = datum0[p].tail(2).normalized();
    datum1[p].tail(2) = datum1[p].tail(2).normalized();
    datum2[p].tail(2) = datum2[p].tail(2).normalized();
  }
  gt_cam[1].block<3,3>(0,0) *= gt_cam[0].block<3,3>(0,0).inverse();
  gt_cam[1].block<3,1>(0,3) -= gt_cam[1].block<3,3>(0,0) * gt_cam[0].block<3,1>(0,3);
  gt_cam[2].block<3,3>(0,0) *= gt_cam[0].block<3,3>(0,0).inverse();
  gt_cam[2].block<3,1>(0,3) -= gt_cam[2].block<3,3>(0,0) * gt_cam[0].block<3,1>(0,3);
  gt_cam[0] = Mat34::Identity();
  for (unsigned v = 0; v < 3; v++)
    OPENMVG_LOG_INFO << "\n "<< gt_cam[v];

  for(unsigned p = 0; p < npoints; ++p) {
    if (!trifocal::NormalizedSquaredPointReprojectionOntoOneViewError::Check(gt_cam, datum0[p], datum1[p], datum2[p]))
      return false;
  }
  return true;
}
#endif
#if 1
// Tests trifocal point-error reprojection tangent-error are very low and that chirality
// pass on perfect synthetic data
//
TEST(SEQUENTIAL_SFM, Trifocal_Check)
{
  NViewOrientedDataSet d;
  nViewDatasetConfigurator config;
  NOrientedPointsCamerasSphere(4, 5, &d, &config); // need to use config since K_ is
                                                    // ignored by getInputScene below
                                                    // (this comes from code
                                                    // initializing K not from a
                                                    // matrix)
  {
  constexpr int ci[3] = {0, 1, 2}; // camera index of the trifocal triplet. Here we
                                   // are skiping camera 0 and using the last 3.
  EXPECT_TRUE(check_camera_triplet(d, ci));
  }
  {
  constexpr int ci[3] = {1, 0, 2}; // camera index of the trifocal triplet. Here we
                                   // are skiping camera 0 and using the last 3.
  EXPECT_TRUE(check_camera_triplet(d, ci));
  }
  {
  constexpr int ci[3] = {2, 0, 1}; // camera index of the trifocal triplet. Here we
                                   // are skiping camera 0 and using the last 3.
  EXPECT_TRUE(check_camera_triplet(d, ci));
  }
  {
  constexpr int ci[3] = {2, 1, 0}; // camera index of the trifocal triplet. Here we
                                   // are skiping camera 0 and using the last 3.
  EXPECT_TRUE(check_camera_triplet(d, ci));
  }
  {
  constexpr int ci[3] = {0, 1, 3}; // camera index of the trifocal triplet. Here we
                                   // are skiping camera 0 and using the last 3.
  EXPECT_TRUE(check_camera_triplet(d, ci));
  }
  {
  constexpr int ci[3] = {0, 2, 3};
  EXPECT_TRUE(check_camera_triplet(d, ci));
  }
  {
  constexpr int ci[3] = {1, 2, 3};
  EXPECT_TRUE(check_camera_triplet(d, ci));
  }
}
#endif


#if 0
// Test a scene where all the camera intrinsics are known
// and oriented features are used for SfM
TEST(SEQUENTIAL_SFM, OrientedSfM) 
{
  const int nviews = synth_nviews_;
  const int npoints = synth_npts_;
  nViewDatasetConfigurator config;
  NViewOrientedDataSet d;
  NOrientedPointsCamerasSphere(nviews, npoints, &d, &config); // need to use config since K_ is
                                              // ignored by getInputScene below
                                              // (this comes from code
                                              // initializing K not from a
                                              // matrix)

  // Translate the input dataset to a SfM_Data scene
  const SfM_Data sfm_data = getInputScene(d, config, PINHOLE_CAMERA);

  // Remove poses and structure
  SfM_Data sfm_data_2 (sfm_data);
  sfm_data_2.poses.clear();
  sfm_data_2.structure.clear();

  SequentialSfMReconstructionEngine sfmEngine(
    sfm_data_2,
   "./",
    stlplus::create_filespec("./", "Reconstruction_Report.html"));

  // Configure the features_provider & the matches_provider from the synthetic dataset
  std::shared_ptr<Features_Provider> feats_provider =
    std::make_shared<Synthetic_Oriented_Features_Provider>();
  // Add a tiny noise in 2D observations to make data more realistic
  dynamic_cast<Synthetic_Oriented_Features_Provider*>(feats_provider.get())->load(d);

  std::shared_ptr<Matches_Provider> matches_provider =
    std::make_shared<Synthetic_Matches_Provider>();
  dynamic_cast<Synthetic_Matches_Provider*>(matches_provider.get())->load(d);
  // Configure data provider (Features and Matches)
  sfmEngine.SetFeaturesProvider(feats_provider.get());
  sfmEngine.SetMatchesProvider(matches_provider.get());
//  XXX sfmEngine.SetResectionMethod(static_cast<resection::SolverType>(static_cast<int>(resection::SolverType::P2Pt_FABBRI_ECCV12)));
  // Configure reconstruction parameters (intrinsic parameters are held constant)
  sfmEngine.Set_Intrinsics_Refinement_Type(cameras::Intrinsic_Parameter_Type::NONE);

  // Will use view ids (1,2,3) as the initial triplet, not (0,1,2)
  assert(nviews > 3); // assuming 4 views
  sfmEngine.setInitialTriplet({sfm_data_2.GetViews().at(1)->id_view,
                               sfm_data_2.GetViews().at(2)->id_view,
                               sfm_data_2.GetViews().at(3)->id_view});
  sfmEngine.SetMaximumTrifocalRansacIterations(1);
  EXPECT_TRUE(sfmEngine.Process());

  const double dResidual = RMSE(sfmEngine.Get_SfM_Data());
  std::cout << "RMSE residual: " << dResidual << std::endl;
  EXPECT_TRUE(dResidual < 0.5);
  EXPECT_TRUE(sfmEngine.Get_SfM_Data().GetPoses().size() == nviews);
  EXPECT_TRUE(sfmEngine.Get_SfM_Data().GetLandmarks().size() == npoints);
  EXPECT_TRUE(IsTracksOneCC(sfmEngine.Get_SfM_Data()));
}
#endif

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
