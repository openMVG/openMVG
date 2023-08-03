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
  template <typename NoiseGenerator>
  bool load( const NViewOrientedDataSet & synthetic_data) {
    // For each view
    for (size_t v = 0; v < synthetic_data._n; ++v) {
      // For each new point visibility
      for (Mat2X::Index i = 0; i < synthetic_data._x[v].cols(); ++i) {
        const Vec2 pt = synthetic_data._x[v].col(i);
        // TODO: fill sio_feats_per_view XXX
        feats_per_view[v].emplace_back(pt(0), pt(1));
      }
    }
    return true;
  }
};


const double K_[2][3] = {
  {2584.9325098195013197, 0, 249.77137587221417903},
  {0, 2584.7918606057692159, 278.31267937919352562}
 //  0 0 1 
};

static constexpr unsigned synth_nviews_ = 4;
static constexpr unsigned synth_npts_ = 4;

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
  { // camera for frame 42
    {-0.097305153950172085242, -0.22322794404612877894, -0.96989741313794208821},
    {0.96072075769186959793, 0.23341709945525662695, -0.15010690664274928263},
    {0.25989869710080021337, -0.94640675329312473618, 0.1917470327448986267},
    {-295.2090359311167731, 1074.0075457376335635, -236.40439390871563319}
  },
  { // camera for frame 54
    {0.81820480546085894158, 0.07511341824191355987, -0.56999900940332604016},
    {0.54313649506229122466, -0.42609616539484057585, 0.72349485524588375007},
    {-0.18853022052765816552, -0.90155423144469115648, -0.38943076883056443327},
    {194.82952402681169701, 1020.3676638972305, 431.76461692675769655}
  },
  { // camera for frame 62
    {-0.61853492444688140672, -0.60388598633423984374, -0.50272881631015808868},
    {-0.025677143402306701336, -0.62392573537516660132, 0.78106168837247091918},
    {-0.78533765448129877473, 0.4960225723146879373, 0.37041379052099593361},
    {887.07508137499985423, -562.68690102473453862, -415.57529638919055515}
  },
  { // camera for frame 07
    {0.42835179674477519285, 0.85688223038479038873, 0.28682325825551502341},
    {0.90183082694967608983, -0.42531965388354808777, -0.076186295248175328609},
    {0.056708886329976546103, 0.29130059263785545998, -0.95494924835828209897},
    {-49.393344597510420613, -315.93153825622943032, 1084.1516821162344968}
  }
};

// Input points and tangents extracted from the synthcurves dataset:
//
// http://github.com/rfabbri/synthcurves-multiview-3d-dataset
//
// Dataset:
// synthcurves-multiview-3d-dataset/spherical-ascii-100_views-perturb-radius_sigma10-normal_sigma0_01rad-minsep_15deg-no_two_cams_colinear_with_object
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
{286.7673976130331539, 217.06531260627261304},
{101.39666157884776965, 215.14917056757076352},
{141.72875899475818073, 270.26093945765563831},
{239.89822517853363593, 86.442049763307068133}
},
// 2D points for frame 54
{
{257.04360648826406077, 159.4404341695463927},
{169.94505787568755295, 309.0591357346235668},
{241.57819514976912956, 446.47561935119495047},
{123.95973916849976604, 213.90676875312345828}
},
// 2D points for frame 62
{
{295.57132984990698787, 147.80261937455236421},
{192.09496151614652604, 285.53721950799990736},
{241.4111595966040511, 409.61663435202348182},
{375.60750199363729962, 277.22372936832925916}
},
// 2D points for frame 07
{
{220.3213713953359445, 152.62347873190785208},
{283.3757713837038068, 153.37782230592358701},
{156.33050901886937822, 334.59959581469149725},
{107.34921953605439171, 97.931489239201212627}
}
};

// 2D tangents for frame 54
// 2D tangents for frame 62
const double 
t_gt_[synth_nviews_][synth_npts_][2] = {
// 2D tangents for frame 42
{
{-0.98898712605989902436, -0.14800156920715851205},
{0.084889345201871871383, -0.99639038487492304075},
{0.9763621177972874321, -0.21614119211847546143},
{0.0082601187924503903515, -0.99996588463683822035}
},
// 2D tangents for frame 54
{
{-0.62769494613031084906, 0.77845941101798377115},
{-0.8302720026964797162, -0.55735841389394735756},
{0.28022966539720772783, -0.95993298444806518521},
{-0.99542450475950383648, 0.095551322985590561587}
},
// 2D tangents for frame 62
{
{-0.53715961235191322931, 0.84348061676480867721},
{0.99931929301637889562, 0.036891064029715384121},
{0.75353934663171129316, -0.65740280884542434681},
{0.76492323888624347283, 0.64412144709812224619}
},
// 2D tangents for frame 07
{
{0.91384373629299286979, -0.406066036055790891},
{-0.42739764583238498696, -0.90406374351421436852},
{-0.91490800607508482312, -0.40366240897526328713},
{-0.89826425202421977811, -0.43945572420366652011}
}
};

const double pts3d_gt_[synth_npts_][3] = {
{-39.999960000000001514, 40.000016999999999712, -39.999979999999993652},
{-28.799995999999886465, 40.000010000000003174, 40.000010000000003174},
{16.241229516856364512, -45.185185185185176238, 41.368080573302677294},
{-83.024179089510298013, -7.2456979436932478222, -4.412526863075626693}
};

const double t3d_gt_[synth_npts_][3] = {
// 3D tangents
{0, 0, 1},
{-1, 0, 0},
{-0.34011103186525448727, -0.10551104075352309153, -0.93444738015720296698},
{-0.71448475613149997621, -0.66437364477423688225, 0.21936087478196092393}
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


/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
