// Copyright (c) 2019 Pierre MOULON.
//:\file
//\Trifocal test with hardcoded samples
//\author Pierre MOULON
//\author Gabriel ANDRADE Rio de Janeiro State U.
//\author Ricardo Fabbri, Brown & Rio de Janeiro State U. (rfabbri.github.io) 

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//
//------------------------------------------------------------------------------
//  This is mostly in openmvg/*/multiview/trifocal/*test*
//  What remains are only sketches of possible tests eg ACRANSAC
//------------------------------------------------------------------------------

#include <testing/testing.h>
#include "openMVG/robust_estimation/robust_estimator_ACRansac.hpp"


#if 0
#define Float double
typedef std::complex<Float> complex;

using namespace trifocal3pt;
using trifocal_model_t = Trifocal3PointPositionTangentialSolver::trifocal_model_t;
typedef MiNuS::minus_util<Float> util;
// using namespace MiNuS;

trifocal_model_t tt_gt_; // corresp. to minus' cameras_gt_

// Converts a trifocal_model to quaternion-translation format
// Assumes tt[0] is identity
void
tt2qt(const trifocal_model_t &tt, double tt_qt[M::nve])
{
  util::rotm2quat(tt[1].transpose().data(), tt_qt);
  util::rotm2quat(tt[2].transpose().data(), tt_qt+4);
  for (unsigned i=0; i < 3; ++i) {
    tt_qt[8+i]   = tt[1](i,3);
    tt_qt[8+3+i] = tt[2](i,3);
  }
}

// Finds a ground truth camera among a std::vector of possible ones
static bool
probe_solutions(
    const std::vector<trifocal_model_t> &solutions, 
    trifocal_model_t &gt, 
    unsigned *solution_index)
{
  Float cameras_quat[M::nsols][M::nve];

  // translate trifocal_model from RC to QT (quaternion-translation)
  // - for each solution
  // -   translate to internal matrix form
  // -   call RC to QT

  tt2qt(gt, data::cameras_gt_quat_);
  for (unsigned s=0; s < solutions.size(); ++s)
    tt2qt(solutions[s], cameras_quat[s]);
  
  return io14::probe_all_solutions_quat(cameras_quat, data::cameras_gt_quat_, solutions.size(), solution_index);
}

static void
initialize_gt()
{
  data::initialize_gt();
  
  double cameras_gt_relative[2][4][3];
  // get relative cameras in usual format
  io14::solution2cams(data::cameras_gt_quat_, cameras_gt_relative);

  tt_gt_[0] = Mat34::Identity(); // view 0 [I | 0]
  for (unsigned v=1; v < io::pp::nviews; ++v) {
    // eigen is col-major but minus is row-major, so memcpy canot be used.
    for (unsigned ir = 0; ir < 3; ++ir)
      for (unsigned ic = 0; ic < 3; ++ic)
        tt_gt_[v](ir, ic) = cameras_gt_relative[v-1][ir][ic];
    for (unsigned r=0; r < 3; ++r) // copy translation
      tt_gt_[v](r,3) = cameras_gt_relative[v-1][3][r];
  }
}
  

TEST(TrifocalSampleApp, solveACRansac) 
{
  {
  std::cerr << "----------------------------------------------------------------\n";
  std::cerr << "3 perfect points = 3 inliers\n";
  array<Mat, 3> datum;   // x,y,orientation across 3 views in normalized world units
  array<Mat, 3> pxdatum; // x,y,orientation across 3 views in pixel units
                         // datum[view](coord,point)
  
  // todo: invert K matrix
  for (unsigned v=0; v < io::pp::nviews; ++v) {
    datum[v].resize(4, 3);
    pxdatum[v].resize(4, 3);
    for (unsigned ip=0; ip < io::pp::npoints; ++ip) {
      datum[v](0,ip) = data::p_[v][ip][0];
      datum[v](1,ip) = data::p_[v][ip][1];
      datum[v](2,ip) = data::tgt_[v][ip][0];
      datum[v](3,ip) = data::tgt_[v][ip][1];
      pxdatum[v] = datum[v];
      trifocal3pt::invert_intrinsics(data::K_, pxdatum[v].col(ip).data(), datum[v].col(ip).data()); 
      trifocal3pt::invert_intrinsics_tgt(data::K_, pxdatum[v].col(ip).data()+2, datum[v].col(ip).data()+2);
    }
  }
  using TrifocalKernel = ThreeViewKernel<Trifocal3PointPositionTangentialSolver, 
                         Trifocal3PointPositionTangentialSolver>;
  
  const TrifocalKernel trifocal_kernel(datum[0], datum[1], datum[2], pxdatum[0], pxdatum[1], pxdatum[2], data::K_);
  
  double threshold = threshold_pixel_to_normalized(1.0, data::K_);
  threshold *= threshold; // squared error
  unsigned constexpr max_iteration = 2; // testing
  // Vector of inliers for the best fit found
  vector<uint32_t> vec_inliers;


  
//  const auto model = MaxConsensus(trifocal_kernel, 
//      ScorerEvaluator<TrifocalKernel>(threshold), &vec_inliers, max_iteration);

  trifocal_model_t model;
  const auto ac_ransac_output = robust::ACRANSAC(
    trifocal_kernel, vec_inliers,
    max_iteration, nullptr,
    std::numeric_limits<double>::infinity(), false);
  
  std::cerr << "Number of inliers (expect 3): "  << vec_inliers.size() << "\n";
  CHECK(vec_inliers.size() == 3);
  }
}
#endif

int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
