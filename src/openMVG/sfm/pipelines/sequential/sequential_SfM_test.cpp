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

// A N-view metric dataset with feature orientation in 3D and 2D.
// All points are seen by all cameras.
struct NViewOrientedDataSet : public NViewDataSet {
  Mat3X _Tgt3d;          // 3D tangent orientation as unit 3D vector.
  std::vector<Mat2X> _tgt2d;  // Projected tangents as unit 2D vector
};

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


const double K_[2][3] = {
  {2584.9325098195013197, 0, 249.77137587221417903},
  {0, 2584.9325098195013197, 278.31267937919352562}
  //  0 0 1 
};

static constexpr unsigned synth_nviews_ = 4;
static constexpr unsigned synth_npts_ = 5;

// camera format: just like a 3x4 [R|T] but transposed to better fit row-major:
// | R |
// | - |
// | T'|
//
// In this case, this data is from the synthcurves multiview dataset,
// so that instead of T, C is stored:
// | R |
// | - |
// | C'|
const double
cameras_gt_[synth_nviews_][4][3] = {
// extrinsics for frame 42
{
{0.032488343069021832776, 0.14118885304658673752, 0.98944945062394118462},
{-0.28679990702465507635, 0.94965585180689460199, -0.12609352267095536027},
{-0.95743946069465868387, -0.27967744082122686367, 0.071345695037685272211},
{1071.484198049582119, 320.76850873549409471, -85.986368935484179588}
},
// extrinsics for frame 54
{
{0.47790123270283219048, -0.019299093028001035321, -0.87820154679288175981},
{-0.29465811799137531235, -0.94535498497159320408, -0.13957272617219795841},
{-0.82751858304384573461, 0.32547119288449349872, -0.45747294709026314896},
{925.05253488236451176, -384.00318581806010343, 531.87270782843597772}
},
// extrinsics for frame 62
{
{-0.219610930168113061, -0.33532970351451174551, -0.91614683828061438398},
{-0.50550514340331198504, -0.76406316361831605466, 0.40083915975658834796},
{-0.83440732819378671259, 0.55114559958548192675, -0.0017142677930838123856},
{951.80671923514557875, -619.60363267357240602, 4.2313312789905133116}
},
// extrinsics for frame 07
{
{0.91074869806248703874, 0.30990046893594913602, -0.27294414873882894002},
{0.29070436143147293517, -0.95055540259185944407, -0.10924926017208461126},
{-0.29330493214776437449, 0.020152567000437736355, -0.95580651327613819213},
{342.08616590607340413, -28.455793406432913883, 1067.8273738052921544}
}
};

// Input points and tangents extracted from the synthcurves dataset:
//
// http://github.com/rfabbri/synthcurves-multiview-3d-dataset
//
// sub-dataset:
// spherical-ascii-100_views-perturb-radius_sigma10-normal_sigma0_01rad-minsep_15deg-no_two_cams_colinear_with_object-aspect_ratio1
//      
//
// Frame files: 42, 54, 62, 07
//
// Points and tangents:
// 
// 620
// 1009
// 3011
// 3389
// 0-based ids. +1 to get file line
// 
// Extracting this data from those files: use scripts/getlines from minus
// or just by hand
//
// This is in pixel image coordinates

const double
p_gt_[synth_nviews_][synth_npts_][2] = {
// 2D points for frame 42
{
{181.53712861181011817, 382.84265753788542952},
{361.39404494701216208, 353.17104859076766843},
{342.08123422244102585, 137.87448982117351193},
{243.09524202092040923, 298.35373638828008325},
{285.88966212823157775, 251.48973104783391364}
},
// 2D points for frame 54
{
{320.61083811349851658, 199.23585641086629039},
{177.83962742245475397, 163.63860158201131867},
{225.74947198803425863, 316.24347634112569949},
{210.92593414503781446, 312.1127002295278885},
{247.68819285385683315, 263.17278766241219046}
},
// 2D points for frame 62
{
{330.38724148135435144, 234.16270784506826885},
{165.33058499047140799, 291.56955369118014687},
{199.86675126597054941, 393.58510880586692338},
{313.99181820108196916, 389.73404770358069982},
{248.50928647922813752, 333.51852292954941959}
},
// 2D points for frame 07
{
{195.98098993753490049, 156.78341685173867859},
{164.52669179292479384, 134.9377776538624687},
{197.46123952507306853, 358.3286736009371225},
{52.110592924391973213, 218.51730104909975694},
{177.00162990735412905, 256.9990394819891435}
}
};

// 2D tangents for frame 54
// 2D tangents for frame 62
const double 
t_gt_[synth_nviews_][synth_npts_][2] = {
// 2D tangents for frame 42
{
{0.99164213923671906681, -0.12901886563609127334},
{-0.27407482149753986667, 0.96170837171207557148},
{-0.99085473006575885968, 0.1349329607853933799},
{0.21688457417910370073, -0.97619725541672619507},
{-0.88923826058271226991, 0.45744433094731007383}
},
// 2D tangents for frame 54
{
{-0.98462682094023079582, -0.1746711867628284176},
{-0.80826947668316584394, 0.58881274872604572046},
{0.90119839239834154121, 0.43340680375213874731},
{-0.54078671415378842813, 0.84115975283815680452},
{0.99935984779679032375, 0.035775614761675351982}
},
// 2D tangents for frame 62
{
{-0.91615034676027939931, 0.4008348065363335766},
{0.44184360558182206313, 0.89709209572175763192},
{0.99085696151653412933, -0.13491657353424593713},
{0.17976067146968616184, 0.98371037454769549857},
{0.91927766178284631149, -0.39360968045395255954}
},
// 2D tangents for frame 07
{
{-0.88483975537184744731, -0.46589548969000454948},
{-0.95661607022265571221, -0.29135149594907394643},
{-0.67890121232474576196, 0.73422962614157050165},
{-0.91684145433395725089, 0.39925148417356498554},
{-0.017735947168850744321, -0.99984270571826627805}
}
};

const double pts3d_gt_[synth_npts_][3] = {
{-39.999960000000001514, 40.000016999999999712, -39.999979999999993652},
{-28.799995999999886465, 40.000010000000003174, 40.000010000000003174},
{16.241229516856364512, -45.185185185185176238, 41.368080573302677294},
{-83.024179089510298013, -7.2456979436932478222, -4.412526863075626693},
{-15.900733289698191442, -6.9926202045388530237, 12.583214033874593696}
};

const double t3d_gt_[synth_npts_][3] = {
{0, 0, 1},
{-1, 0, 0},
{-0.34011103186525448727, -0.10551104075352309153, -0.93444738015720296698},
{-0.71448475613149997621, -0.66437364477423688225, 0.21936087478196092393},
{-0.38185692861567399614, 0.23333310127926898403, -0.89428236587534382096}
};

// number of points is hardcoded and number of views is hardcoded
void 
NOrientedPointsCamerasSphere(NViewOrientedDataSet *dp)
{
  //-- Setup a camera rig
  // based on github.com/rfabbri/synthcurves-multiview-3d-dataset
  // by hardcoding points
  NViewOrientedDataSet &d = *dp;

  unsigned const nviews = synth_nviews_;
  unsigned const npoints = synth_npts_;
  d._n = nviews;
  d._K.resize(nviews);
  d._R.resize(nviews);
  d._t.resize(nviews);
  d._C.resize(nviews);
  d._x.resize(nviews);
  d._tgt2d.resize(nviews);
  d._x_ids.resize(nviews);
  d._X.resize(3, npoints);
  d._Tgt3d.resize(3, npoints);

  Vecu all_point_ids(npoints);
  for (size_t p = 0; p < npoints; ++p) {
    all_point_ids[p] = p;
    d._X.col(p) << pts3d_gt_[p][0], pts3d_gt_[p][1], pts3d_gt_[p][2];
    d._Tgt3d.col(p) << t3d_gt_[p][0], t3d_gt_[p][1], t3d_gt_[p][2];
  }

  for (size_t v = 0; v < nviews; ++v) {
    d._C[v] << cameras_gt_[v][3][0] , cameras_gt_[v][3][1] , cameras_gt_[v][3][2];

    d._K[v] << K_[0][0],           K_[0][1], K_[0][2],
               K_[1][0],           K_[1][1], K_[1][2],
                      0,           0,          1;
    d._R[v] << cameras_gt_[v][0][0] , cameras_gt_[v][0][1] , cameras_gt_[v][0][2],
               cameras_gt_[v][1][0] , cameras_gt_[v][1][1] , cameras_gt_[v][1][2],
               cameras_gt_[v][2][0] , cameras_gt_[v][2][1] , cameras_gt_[v][2][2];

    d._t[v] = -d._R[v] * d._C[v];
    d._x[v].resize(2,npoints);
    d._tgt2d[v].resize(2,npoints);
    for (unsigned p = 0; p < npoints; ++p) {
      d._x[v].col(p) << p_gt_[v][p][0], p_gt_[v][p][1];
      d._tgt2d[v].col(p) << t_gt_[v][p][0], t_gt_[v][p][1];
    }
    d._x_ids[v] = all_point_ids;
  }
}


// Test a scene where all the camera intrinsics are known
// and oriented features are used for SfM
TEST(SEQUENTIAL_SFM, OrientedSfM) {

  const int nviews = synth_nviews_;
  const int npoints = synth_npts_;
  const nViewDatasetConfigurator config;
  NViewOrientedDataSet d;
  NOrientedPointsCamerasSphere(&d);

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
    std::make_shared<Synthetic_Oriented_Features_Provider>();
  // Add a tiny noise in 2D observations to make data more realistic
  dynamic_cast<Synthetic_Oriented_Features_Provider*>(feats_provider.get())->load(d);

  std::shared_ptr<Matches_Provider> matches_provider =
    std::make_shared<Synthetic_Matches_Provider>();
  dynamic_cast<Synthetic_Matches_Provider*>(matches_provider.get())->load(d);

  // Configure data provider (Features and Matches)
  sfmEngine.SetFeaturesProvider(feats_provider.get());
  sfmEngine.SetMatchesProvider(matches_provider.get());

  // Configure reconstruction parameters (intrinsic parameters are held constant)
  sfmEngine.Set_Intrinsics_Refinement_Type(cameras::Intrinsic_Parameter_Type::NONE);

  // Will use view ids (1,2,3) as the initial triplet
  assert(nviews > 3); // assuming 4 views
  sfmEngine.setInitialTriplet({sfm_data_2.GetViews().at(1)->id_view,
                            sfm_data_2.GetViews().at(2)->id_view,
                            sfm_data_2.GetViews().at(3)->id_view});
  sfmEngine.SetMaximumTrifocalRansacIterations(5);
  EXPECT_TRUE (sfmEngine.Process());

  const double dResidual = RMSE(sfmEngine.Get_SfM_Data());
  std::cout << "RMSE residual: " << dResidual << std::endl;
  EXPECT_TRUE( dResidual < 0.5);
  EXPECT_TRUE( sfmEngine.Get_SfM_Data().GetPoses().size() == nviews);
  EXPECT_TRUE( sfmEngine.Get_SfM_Data().GetLandmarks().size() == npoints);
  EXPECT_TRUE( IsTracksOneCC(sfmEngine.Get_SfM_Data()));
}


/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
