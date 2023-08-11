// Copyright (c) 2019 Pierre MOULON.
//:\file
//\Trifocal test with hardcoded samples
//\author Pierre MOULON
//\author Gabriel ANDRADE Rio de Janeiro State U.
//\author Ricardo Fabbri, Brown & Rio de Janeiro State U. (rfabbri.github.io) 

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <testing/testing.h>
#include <minus/chicago-default.h>
#include "openMVG/multiview/trifocal/solver_trifocal_metrics.hpp"

using namespace std;
using namespace openMVG;
using namespace openMVG::trifocal;

trifocal_model_t tt_gt_; // corresp. to minus' cameras_gt_

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

static void
invert_intrinsics_tgt(
    const double K[/*3 or 2 ignoring last line*/][3], 
    const double px_tgt_coords[2], 
    double normalized_tgt_coords[2])
{
  const double *tp = px_tgt_coords;
  double *t = normalized_tgt_coords;
  t[1] = tp[1]/K[1][1];
  t[0] = (tp[0] - K[0][1]*t[1])/K[0][0];
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

// Testing Error()
TEST(TrifocalSampleApp, error_simple) 
{
  initialize_gt(); // gets tt_gt_
  
  {
  std::cerr << "----------------------------------------------------------------\n";
  std::cerr << "Error model with 3 perfect points\n";
  
  array<Mat, 3> datum; // x,y,orientation across 3 views
  array<Mat, 3> pxdatum; // for debug
  
  for (unsigned v=0; v < 3; ++v) {
    datum[v].resize(4, 1);
    pxdatum[v].resize(4, 1);
    constexpr unsigned ip = 0;
    pxdatum[v](0,ip) = data::p_[v][ip][0];
    pxdatum[v](1,ip) = data::p_[v][ip][1];
    pxdatum[v](2,ip) = data::tgt_[v][ip][0];
    pxdatum[v](3,ip) = data::tgt_[v][ip][1];
    invert_intrinsics(data::K_, pxdatum[v].col(ip).data(), datum[v].col(ip).data()); 
    invert_intrinsics_tgt(data::K_, pxdatum[v].col(ip).data()+2, datum[v].col(ip).data()+2);
  }
  
  float err = NormalizedSquaredPointReprojectionOntoOneViewError::Error(tt_gt_, 
      datum[0].col(0), datum[1].col(0), datum[2].col(0)); 
  
  std::cerr << "Error (squared, normalized): " << err << "\n";
  std::cerr << "Error (pixel, not squared): " << 
    NormalizedSquaredPointReprojectionOntoOneViewError::threshold_normalized_to_pixel(sqrt(err),data::K_) << "\n";
  CHECK(err < 1e-6);

  CHECK(NormalizedSquaredPointReprojectionOntoOneViewError::Check(tt_gt_, datum[0].col(0), datum[1].col(0), datum[2].col(0)));
  }
  
  // Testing error model with 2 perfect points and 1 perturbed point
  {
  std::cerr << "----------------------------------------------------------------\n";
  std::cerr << "Testing error model with 2 perfect points and 1 perturbed point\n";
  array<Mat, 3> datum; // x,y,orientation across 3 views
  array<Mat, 3> pxdatum;
  
  for (unsigned v = 0; v < 3; ++v) {
    datum[v].resize(4,1);
    pxdatum[v].resize(4,1);
    constexpr unsigned ip = 0;
    pxdatum[v](0,ip) = data::p_[v][ip][0];
    pxdatum[v](1,ip) = data::p_[v][ip][1];
    pxdatum[v](2,ip) = data::tgt_[v][ip][0];
    pxdatum[v](3,ip) = data::tgt_[v][ip][1];
    invert_intrinsics(data::K_, pxdatum[v].col(ip).data(), datum[v].col(ip).data()); 
    invert_intrinsics_tgt(data::K_, pxdatum[v].col(ip).data()+2, datum[v].col(ip).data()+2);
  }

  {
    constexpr unsigned v = 1, ip = 0;
    pxdatum[v](0,ip) += 5;
    pxdatum[v](1,ip) += 5;
    invert_intrinsics(data::K_, pxdatum[v].col(ip).data(), datum[v].col(ip).data()); 
  }

  float err = NormalizedSquaredPointReprojectionOntoOneViewError::Error(tt_gt_, 
      datum[0].col(0), datum[1].col(0), datum[2].col(0)); 
  
  std::cerr << "Error (squared, normalized): " << err << "\n";
  std::cerr << "Error (pixel, not squared): " << 
    NormalizedSquaredPointReprojectionOntoOneViewError::threshold_normalized_to_pixel(sqrt(err),data::K_) << "\n";
  CHECK(NormalizedSquaredPointReprojectionOntoOneViewError::threshold_normalized_to_pixel(sqrt(err),data::K_) > 1);
  }
}

int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
