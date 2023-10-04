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

#include "openMVG/robust_estimation/robust_estimator_MaxConsensus.hpp"
#include "openMVG/robust_estimation/score_evaluator.hpp"
#include "openMVG/multiview/trifocal/solver_trifocal_three_point.hpp"
#include "openMVG/multiview/trifocal/solver_trifocal_metrics.hpp"
#include "openMVG/multiview/trifocal/three_view_kernel.hpp"

using namespace std;
using namespace openMVG;
using namespace openMVG::trifocal;
using namespace openMVG::robust;

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
  double n = hypot(t[0],t[1]);
  t[0] /= n; t[1] /= n;
}

TEST(TrifocalSampleApp, solveRansac)
{
  {
  std::cerr << "----------------------------------------------------------------\n";
  std::cerr << "3 perfect points = 3 inliers\n";
  array<Mat, 3> datum;   // x,y,orientation across 3 views in normalized world units

  for (unsigned v=0; v < io::pp::nviews; ++v) {
    datum[v].resize(4, 3);
    for (unsigned ip=0; ip < io::pp::npoints; ++ip) {
      datum[v](0,ip) = data::p_[v][ip][0];
      datum[v](1,ip) = data::p_[v][ip][1];
      datum[v](2,ip) = data::tgt_[v][ip][0];
      datum[v](3,ip) = data::tgt_[v][ip][1];
      invert_intrinsics(data::K_, datum[v].col(ip).data(), datum[v].col(ip).data());
      invert_intrinsics_tgt(data::K_, datum[v].col(ip).data()+2, datum[v].col(ip).data()+2);
    }
  }
  using TrifocalKernel = ThreeViewKernel<Trifocal3PointPositionTangentialSolver,
                         NormalizedSquaredPointReprojectionOntoOneViewError>;

  const TrifocalKernel trifocal_kernel(datum[0], datum[1], datum[2]);

  double threshold =
    NormalizedSquaredPointReprojectionOntoOneViewError::threshold_pixel_to_normalized(1e-5, data::K_);
  threshold *= threshold; // squared error
  unsigned constexpr max_iteration = 2; // testing
  // Vector of inliers for the best fit found
  vector<uint32_t> vec_inliers;
  const auto model = MaxConsensus(trifocal_kernel,
      ScorerEvaluator<TrifocalKernel>(threshold), &vec_inliers, max_iteration);
  std::cerr << "Number of inliers (expect 3): "  << vec_inliers.size() << "\n";
  CHECK(vec_inliers.size() == 3);
  }

  {
  std::cerr << "----------------------------------------------------------------\n";
  std::cerr << "3 perfect points and 1 outlier\n";
  array<Mat, 3> datum;   // x,y,orientation across 3 views in normalized world units

  for (unsigned v=0; v < io::pp::nviews; ++v) {
    datum[v].resize(4, 4);
    for (unsigned ip=0; ip < io::pp::npoints; ++ip) {
      datum[v](0,ip) = data::p_[v][ip][0];
      datum[v](1,ip) = data::p_[v][ip][1];
      datum[v](2,ip) = data::tgt_[v][ip][0];
      datum[v](3,ip) = data::tgt_[v][ip][1];
      invert_intrinsics(data::K_, datum[v].col(ip).data(), datum[v].col(ip).data());
      invert_intrinsics_tgt(data::K_, datum[v].col(ip).data()+2, datum[v].col(ip).data()+2);
    }
    // 4th point is just the 2nd one, perturbed
    constexpr unsigned ip = 3;
    datum[v](0,ip) = data::p_[v][1][0];
    datum[v](1,ip) = data::p_[v][1][1];
    datum[v](2,ip) = data::tgt_[v][1][0];
    datum[v](3,ip) = data::tgt_[v][1][1];
    datum[v](0,ip) += 5.0;
    datum[v](1,ip) += 5.0;
    invert_intrinsics(data::K_, datum[v].col(ip).data(), datum[v].col(ip).data());
    invert_intrinsics_tgt(data::K_, datum[v].col(ip).data()+2, datum[v].col(ip).data()+2);
  }

  using TrifocalKernel = ThreeViewKernel<Trifocal3PointPositionTangentialSolver,
                         NormalizedSquaredPointReprojectionOntoOneViewError>;

  const TrifocalKernel trifocal_kernel(datum[0], datum[1], datum[2]);

  double threshold =
    NormalizedSquaredPointReprojectionOntoOneViewError::threshold_pixel_to_normalized(1.0, data::K_);
  threshold *= threshold; // squared error
  unsigned constexpr max_iteration = 2; // testing
  // Vector of inliers for the best fit found
  vector<uint32_t> vec_inliers;
  const auto model = MaxConsensus(trifocal_kernel,
      ScorerEvaluator<TrifocalKernel>(threshold), &vec_inliers, max_iteration);
  std::cerr << "Number of inliers (expect never 4): "  << vec_inliers.size() << "\n";
  CHECK(vec_inliers.size() <= 3);
  }

  {
  std::cerr << "----------------------------------------------------------------\n";
  std::cerr << "3 perfect points, 1 outlier, increase threshold\n";
  array<Mat, 3> datum;   // x,y,orientation across 3 views in normalized world units

  for (unsigned v=0; v < io::pp::nviews; ++v) {
    datum[v].resize(4, 4);
    for (unsigned ip=0; ip < io::pp::npoints; ++ip) {
      datum[v](0,ip) = data::p_[v][ip][0];
      datum[v](1,ip) = data::p_[v][ip][1];
      datum[v](2,ip) = data::tgt_[v][ip][0];
      datum[v](3,ip) = data::tgt_[v][ip][1];
      invert_intrinsics(data::K_, datum[v].col(ip).data(), datum[v].col(ip).data());
      invert_intrinsics_tgt(data::K_, datum[v].col(ip).data()+2, datum[v].col(ip).data()+2);
    }
    // 4th point is just the 2nd one, perturbed
    constexpr unsigned ip = 3;
    datum[v](0,ip) = data::p_[v][1][0];
    datum[v](1,ip) = data::p_[v][1][1];
    datum[v](2,ip) = data::tgt_[v][1][0];
    datum[v](3,ip) = data::tgt_[v][1][1];
    datum[v](0,ip) += 5.0;
    datum[v](1,ip) += 5.0;
    invert_intrinsics(data::K_, datum[v].col(ip).data(), datum[v].col(ip).data());
    invert_intrinsics_tgt(data::K_, datum[v].col(ip).data()+2, datum[v].col(ip).data()+2);
  }

  using TrifocalKernel = ThreeViewKernel<Trifocal3PointPositionTangentialSolver,
                         NormalizedSquaredPointReprojectionOntoOneViewError>;

  const TrifocalKernel trifocal_kernel(datum[0], datum[1], datum[2]);

  double threshold =
    NormalizedSquaredPointReprojectionOntoOneViewError::threshold_pixel_to_normalized(20.0, data::K_);
  threshold *= threshold; // squared error
  unsigned constexpr max_iteration = 2; // testing
  // Vector of inliers for the best fit found
  vector<uint32_t> vec_inliers;
  const auto model = MaxConsensus(trifocal_kernel,
      ScorerEvaluator<TrifocalKernel>(threshold), &vec_inliers, max_iteration);
  std::cerr << "Number of inliers (expect 4): "  << vec_inliers.size() << "\n";
  CHECK(vec_inliers.size() == 4);
  }

  {
  std::cerr << "----------------------------------------------------------------\n";
  std::cerr << "3 perfect points, 1 outlier, 2 repeated perfect points\n";
  array<Mat, 3> datum;   // x,y,orientation across 3 views in normalized world units

  for (unsigned v=0; v < io::pp::nviews; ++v) {
    datum[v].resize(4, 6);
    for (unsigned ip=0; ip < io::pp::npoints; ++ip) {
      datum[v](0,ip) = data::p_[v][ip][0];
      datum[v](1,ip) = data::p_[v][ip][1];
      datum[v](2,ip) = data::tgt_[v][ip][0];
      datum[v](3,ip) = data::tgt_[v][ip][1];
      invert_intrinsics(data::K_, datum[v].col(ip).data(), datum[v].col(ip).data());
      invert_intrinsics_tgt(data::K_, datum[v].col(ip).data()+2, datum[v].col(ip).data()+2);
    }
    // 4th point is just the 2nd one, perturbed
    {
    constexpr unsigned ip = 3;
    datum[v](0,ip) = data::p_[v][1][0];
    datum[v](1,ip) = data::p_[v][1][1];
    datum[v](2,ip) = data::tgt_[v][1][0];
    datum[v](3,ip) = data::tgt_[v][1][1];
    datum[v](0,ip) += 5.0;
    datum[v](1,ip) += 5.0;
    invert_intrinsics(data::K_, datum[v].col(ip).data(), datum[v].col(ip).data());
    invert_intrinsics_tgt(data::K_, datum[v].col(ip).data()+2, datum[v].col(ip).data()+2);
    }

    // 5th and 6th are 1st repeated
    datum[v].col(4) = datum[v].col(0);
    datum[v].col(5) = datum[v].col(0);
  }

  using TrifocalKernel = ThreeViewKernel<Trifocal3PointPositionTangentialSolver,
                         NormalizedSquaredPointReprojectionOntoOneViewError>;

  const TrifocalKernel trifocal_kernel(datum[0], datum[1], datum[2]);

  double threshold =
    NormalizedSquaredPointReprojectionOntoOneViewError::threshold_pixel_to_normalized(1.0, data::K_);
  threshold *= threshold; // squared error
  unsigned constexpr max_iteration = 5; // testing
  // Vector of inliers for the best fit found
  vector<uint32_t> vec_inliers;
  const auto model = MaxConsensus(trifocal_kernel,
      ScorerEvaluator<TrifocalKernel>(threshold), &vec_inliers, max_iteration);
  std::cerr << "Number of inliers (expect 5): "  << vec_inliers.size() << "\n";
  CHECK(vec_inliers.size() == 5);
  }
}

int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
