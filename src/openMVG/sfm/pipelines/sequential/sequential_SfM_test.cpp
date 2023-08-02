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

// A N-view metric dataset with feature orientation in 3D and 2D.
// All points are seen by all cameras.
struct NViewOrientedDataSet : public NViewDataSet {
  Mat3X _T;          // 3D tangent orientation as unit 3D vector.
  std::vector<Mat2X> _t;  // Projected tangents as unit 2D vector
};


K_[2][3] = {
  {2584.9325098195013197, 0, 249.77137587221417903},
  {0, 2584.7918606057692159, 278.31267937919352562}
 //  0 0 1 
};




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
static constexpr synth_nviews = 3;
template <typename F>
F minus_data<cleveland14a,F>::
cameras_gt_[synth_nviews][4][3] = {
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
  }
};

// Input points and tangents corresponding to the above gammified homotopy parameters,
// extracted from the original synthcurves dataset:
//
// http://github.com/rfabbri/synthcurves-multiview-3d-dataset
// 
// synthcurves-multiview-3d-dataset/spherical-ascii-100_views-perturb-radius_sigma10-normal_sigma0_01rad-minsep_15deg-no_two_cams_colinear_with_object
//
// Frame files: frame_..42, 54, 62
//
// Points:
// 
// 620
// 3011  tangents
// 3389  tangents
// 0-based ids. +1 to get file line
// 
// NOTE: This is exactly the input case in Fabbri's slides in big-notes/trifocal/trifocal.key
// However, other developers made an off by 1 mistake, and the actual points
// obtained for this data are 219 3010 and 3388 using 0-based ids.
// The original one as described in the above keynote slides are in p_correct_
// 
// Extracting this data from those files: use scripts/getlines.sh from
// synthcurves dataset.
//
// This is in pixel image coordinates

const double
p_gt_[io::pp::nviews][io::pp::npoints][io::ncoords2d] = {
  // points for frame 42
  // + sed -n '3012p;3390p;621p' frame_0042-pts-2d.txt
  {
    {141.01103052308988595, 270.45312297462106699},
    {239.89822517853363593, 86.442049763307068133},
    {286.7673976130331539, 217.06531260627261304}
  },
  // points for frame 54
  // + sed -n '3012p;3390p;621p;' frame_0054-pts-2d.txt
  {
    {241.41513314836856807, 447.15662243793082098},
    {123.95973916849976604, 213.90676875312345828},
    {257.04360648826406077, 159.4404341695463927}
  },
  // points for frame 62
  // + sed -n '3012p;3390p;621p' frame_0062-pts-2d.txt
  {
    {375.60750199363729962, 277.22372936832925916},
    {295.57132984990698787, 147.80261937455236421},
    {240.78946527513195974, 410.13737156824942076}
  }
};

// 2D tangents for frame 54
// 2D tangents for frame 62
const double minus_data<cleveland14a,F>::
t_gt_[io::pp::nviews][io::pp::npoints][io::ncoords2d] = {
  // tangents for frame 42
  // + sed -n '3012p;3390p' frame_0042-tgts-2d.txt
  {
    {0.9536809622336909209, -0.3008199166827579818},
    {0.0082601187924503903515, -0.99996588463683822035}
    // zero (unused by us)
  },
  // tangents for frame 54
  // + sed -n '3012p;3390p' frame_0054-tgts-2d.txt
  {
    {0.18491347256048701331, -0.9827548054655455001},
    {-0.99542450475950383648, 0.095551322985590561587}
    // zero (unused by us)
  },
  // tangents for frame 62
  // + sed -n '3012p;3390p' frame_0062-tgts-2d.txt
  {
    {0.77931350598248894102, -0.62663423094599701724},
    {0.76492323888624347283, 0.64412144709812224619}
    // zero (unused by us)
  }
};

const double pts3d_gt_[io::pp::npoints][io::ncoords3d] = {
}

unsigned point_ids_ = {
};


// number of points is hardcoded and number of views is hardcoded
void 
NOrientedPointsCamerasSphere(NViewOrientedDataSet *dp)
{
  //-- Setup a camera rig
  // based on github.com/rfabbri/synthcurves-multiview-3d-dataset
  // by hardcoding points
  NViewOrientedDataSet &d = *dp;

  unsigned nviews = 4;
  unsigned npoints = 6;
  d._n = nviews;
  d._K.resize(nviews);
  d._R.resize(nviews);
  d._t.resize(nviews);
  d._C.resize(nviews);
  d._x.resize(nviews);
  d._x_ids.resize(nviews);
  d._X.resize(3, npoints);

  Vecu all_point_ids(npoints);
  for (size_t p = 0; p < npoints; ++p) {
    all_point_ids[p] = point_ids[p];
    d._X.col(p) = pts3d_gt_[p][0] << pts3d_gt_[p][1] << pts3d_gt_[p][2];
    d._T.col(p) = t3d_gt_[p][0] << t3d_gt_[p][1] << t3d_gt_[p][2];
  }

  for (size_t v = 0; v < nviews; ++v) {
    d._C[v] << cameras_gt_[v][3][1] , cameras_gt_[v][3][2] , cameras_gt_[v][3][3];

    d._K[v] << K_[0][0],           K_[0][1], K_[0][2],
               K_[1][0],           K_[1][1], K_[1][2],
                        0,           0,          1;
    d._R[v] = Eigen::Map<Matrix<double,3,3,RowMajor> >(cameras_gt_);
    d._t[v] = -d._R[v] * camera_center;
    d._x[v].resize(2,npoints);
    d._t[v].resize(2,npoints);
    for (unsigned p = 0; p < npoints; ++p) {
      d._x[v].col(p) << p_gt_[p][0], p_gt_[p][1];
      d._t[v].col(p) << t_gt_[p][0], t_gt_[p][1];
    }
    d._x_ids[v] = all_point_ids;
  }
}


// Test a scene where all the camera intrinsics are known
// and oriented features are used for SfM
TEST(SEQUENTIAL_SFM, OrientedSfM) {

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


/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
