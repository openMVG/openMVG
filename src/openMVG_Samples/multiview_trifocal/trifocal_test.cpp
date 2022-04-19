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
#include "minus/chicago-default.h"
#include "minus/internal-util.h"
#include <minus/debug-common.h>
// #include "minus/debug-common.h"
#include "trifocal.h"
#include "trifocal-app.h"
#include "trifocal-util.h"

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
  

// Directly runs the solver and test
// - define synthetic data
// - directly passo to solver
//    - use Trifocal3PointPositionTangentialSolver::Solve
// - check if the solver returns any known root
// - might fail 5% of the time
// - also check if Trifocal3PointPositionTangentialSolver::error function returns zero
TEST(TrifocalSampleApp, solver) 
{
  array<Mat, 3> datum; // x,y,orientation across 3 views
                       // datum[view](coord,point)
  
  for (unsigned v=0; v < 3; ++v) {
    datum[v].resize(4, 3);
    for (unsigned ip=0; ip < 3; ++ip) {
      datum[v](0,ip) = data::p_[v][ip][0];
      datum[v](1,ip) = data::p_[v][ip][1];
      datum[v](2,ip) = data::tgt_[v][ip][0];
      datum[v](3,ip) = data::tgt_[v][ip][1];
      trifocal3pt::invert_intrinsics(data::K_, datum[v].col(ip).data(), datum[v].col(ip).data()); 
      trifocal3pt::invert_intrinsics_tgt(data::K_, datum[v].col(ip).data()+2, datum[v].col(ip).data()+2);
    }
  }
  
  initialize_gt();
  unsigned constexpr max_solve_tries = 5;
  bool found = false;

  for (unsigned i = 0; i < max_solve_tries; ++i) {
    std::vector<trifocal_model_t> sols; // std::vector of 3-cam solutions
    
    std::cerr << "Test log: Trying to solve, attempt: " << i+1 << std::endl;
    Trifocal3PointPositionTangentialSolver::Solve(datum[0], datum[1], datum[2], &sols);

    unsigned sol_id = (unsigned)-1;
    found = probe_solutions(sols, tt_gt_, &sol_id);
    if (found) {
      std::cerr << "Found solution at id " << sol_id << std::endl;
      //      for(unsigned j = 0; j < 3; j++)
      //        std::cout << sols[sol_id][j] << "\n" << std::endl;
      break;
    }
    std::cerr << "Test log: Solve failed to find ground truth. Retrying different randomization\n";
  }
  
  CHECK(found);
}

// Testing Error()
TEST(TrifocalSampleApp, error_simple) 
{
  initialize_gt(); // gets tt_gt_
  
  {
  std::cerr << "----------------------------------------------------------------\n";
  std::cerr << "Error model with 3 perfect points\n";
  
  array<Mat, 3> datum; // x,y,orientation across 3 views
  array<Mat, 3> pxdatum;
  
  for (unsigned v=0; v < 3; ++v) {
    datum[v].resize(4, 1);
    pxdatum[v].resize(4, 1);
    constexpr unsigned ip = 0;
    pxdatum[v](0,ip) = data::p_[v][ip][0];
    pxdatum[v](1,ip) = data::p_[v][ip][1];
    pxdatum[v](2,ip) = data::tgt_[v][ip][0];
    pxdatum[v](3,ip) = data::tgt_[v][ip][1];
    trifocal3pt::invert_intrinsics(data::K_, pxdatum[v].col(ip).data(), datum[v].col(ip).data()); 
    trifocal3pt::invert_intrinsics_tgt(data::K_, pxdatum[v].col(ip).data()+2, datum[v].col(ip).data()+2);
  }
  
  float err = Trifocal3PointPositionTangentialSolver::Error(tt_gt_, 
      datum[0].col(0), datum[1].col(0), datum[2].col(0), 
      pxdatum[0].col(0), pxdatum[1].col(0), pxdatum[2].col(0), data::K_); 
  
  std::cerr << "Error (squared, normalized): " << err << "\n";
  std::cerr << "Error (pixel, not squared): " << threshold_normalized_to_pixel(sqrt(err),data::K_) << "\n";
  CHECK(err < 0.001);
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
    trifocal3pt::invert_intrinsics(data::K_, pxdatum[v].col(ip).data(), datum[v].col(ip).data()); 
    trifocal3pt::invert_intrinsics_tgt(data::K_, pxdatum[v].col(ip).data()+2, datum[v].col(ip).data()+2);
  }

  {
    constexpr unsigned v = 1, ip = 0;
    pxdatum[v](0,ip) += 5;
    pxdatum[v](1,ip) += 5;
    trifocal3pt::invert_intrinsics(data::K_, pxdatum[v].col(ip).data(), datum[v].col(ip).data()); 
  }

  float err = Trifocal3PointPositionTangentialSolver::Error(tt_gt_, 
      datum[0].col(0), datum[1].col(0), datum[2].col(0), 
      pxdatum[0].col(0), pxdatum[1].col(0), pxdatum[2].col(0), data::K_); 
  
  std::cerr << "Error (squared, normalized): " << err << "\n";
  std::cerr << "Error (pixel, not squared): " << threshold_normalized_to_pixel(sqrt(err),data::K_) << "\n";
  CHECK(threshold_normalized_to_pixel(sqrt(err),data::K_) > 1);
  }
}

#include "openMVG/robust_estimation/robust_estimator_MaxConsensus.hpp"
#include "openMVG/robust_estimation/score_evaluator.hpp"

using namespace openMVG::robust;

// Runs the solve through ransac 
// - first, synthetic data with three perfect points
// - second, synthetic data with one outlier
TEST(TrifocalSampleApp, solveRansac) 
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
  unsigned constexpr max_iteration = 2; // testing
  // Vector of inliers for the best fit found
  vector<uint32_t> vec_inliers;
  const auto model = MaxConsensus(trifocal_kernel, 
      ScorerEvaluator<TrifocalKernel>(threshold), &vec_inliers, max_iteration);
  std::cerr << "Number of inliers (expect 3): "  << vec_inliers.size() << "\n";
  CHECK(vec_inliers.size() == 3);
      
  // TODO check the cameras
  }

  {
  std::cerr << "----------------------------------------------------------------\n";
  std::cerr << "3 perfect points and 1 outlier\n";
  array<Mat, 3> datum;   // x,y,orientation across 3 views in normalized world units
  array<Mat, 3> pxdatum; // x,y,orientation across 3 views in pixel units
                         // datum[view](coord,point)
  
  // todo: invert K matrix
  for (unsigned v=0; v < io::pp::nviews; ++v) {
    datum[v].resize(4, 4);
    pxdatum[v].resize(4, 4);
    for (unsigned ip=0; ip < io::pp::npoints; ++ip) {
      datum[v](0,ip) = data::p_[v][ip][0];
      datum[v](1,ip) = data::p_[v][ip][1];
      datum[v](2,ip) = data::tgt_[v][ip][0];
      datum[v](3,ip) = data::tgt_[v][ip][1];
      pxdatum[v] = datum[v];
      trifocal3pt::invert_intrinsics(data::K_, pxdatum[v].col(ip).data(), datum[v].col(ip).data()); 
      trifocal3pt::invert_intrinsics_tgt(data::K_, pxdatum[v].col(ip).data()+2, datum[v].col(ip).data()+2);
    }
    // 4th point is just the 2nd one, perturbed
    constexpr unsigned ip = 3;
    pxdatum[v].col(ip) = pxdatum[v].col(1);
    pxdatum[v](0,ip) += 5.0;
    pxdatum[v](1,ip) += 5.0;
    trifocal3pt::invert_intrinsics(data::K_, pxdatum[v].col(ip).data(), datum[v].col(ip).data()); 
    trifocal3pt::invert_intrinsics_tgt(data::K_, pxdatum[v].col(ip).data()+2, datum[v].col(ip).data()+2);
  }
  
  using TrifocalKernel = ThreeViewKernel<Trifocal3PointPositionTangentialSolver, 
                         Trifocal3PointPositionTangentialSolver>;
  
  const TrifocalKernel trifocal_kernel(datum[0], datum[1], datum[2], pxdatum[0], pxdatum[1], pxdatum[2], data::K_);
  
  double threshold = threshold_pixel_to_normalized(1.0, data::K_);
  threshold *= threshold; // squared error
  std::cerr << "Threshold squared in normalized coords "  << threshold << std::endl;
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
  array<Mat, 3> pxdatum; // x,y,orientation across 3 views in pixel units
                         // datum[view](coord,point)
  
  // todo: invert K matrix
  for (unsigned v=0; v < io::pp::nviews; ++v) {
    datum[v].resize(4, 4);
    pxdatum[v].resize(4, 4);
    for (unsigned ip=0; ip < io::pp::npoints; ++ip) {
      datum[v](0,ip) = data::p_[v][ip][0];
      datum[v](1,ip) = data::p_[v][ip][1];
      datum[v](2,ip) = data::tgt_[v][ip][0];
      datum[v](3,ip) = data::tgt_[v][ip][1];
      pxdatum[v] = datum[v];
      trifocal3pt::invert_intrinsics(data::K_, pxdatum[v].col(ip).data(), datum[v].col(ip).data()); 
      trifocal3pt::invert_intrinsics_tgt(data::K_, pxdatum[v].col(ip).data()+2, datum[v].col(ip).data()+2);
    }
    // 4th point is just the 2nd one, perturbed
    constexpr unsigned ip = 3;
    pxdatum[v].col(ip) = pxdatum[v].col(1);
    pxdatum[v](0,ip) += 5.0;
    pxdatum[v](1,ip) += 5.0;
    trifocal3pt::invert_intrinsics(data::K_, pxdatum[v].col(ip).data(), datum[v].col(ip).data()); 
    trifocal3pt::invert_intrinsics_tgt(data::K_, pxdatum[v].col(ip).data()+2, datum[v].col(ip).data()+2);
  }
  
  using TrifocalKernel = ThreeViewKernel<Trifocal3PointPositionTangentialSolver, 
                         Trifocal3PointPositionTangentialSolver>;
  
  const TrifocalKernel trifocal_kernel(datum[0], datum[1], datum[2], pxdatum[0], pxdatum[1], pxdatum[2], data::K_);
  
  double threshold = threshold_pixel_to_normalized(100.0, data::K_);
  threshold *= threshold; // squared error
  std::cerr << "Threshold squared in normalized coords "  << threshold << std::endl;
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
  array<Mat, 3> pxdatum; // x,y,orientation across 3 views in pixel units
                         // datum[view](coord,point)
  
  // todo: invert K matrix
  for (unsigned v=0; v < io::pp::nviews; ++v) {
    datum[v].resize(4, 6);
    pxdatum[v].resize(4, 6);
    for (unsigned ip=0; ip < io::pp::npoints; ++ip) {
      datum[v](0,ip) = data::p_[v][ip][0];
      datum[v](1,ip) = data::p_[v][ip][1];
      datum[v](2,ip) = data::tgt_[v][ip][0];
      datum[v](3,ip) = data::tgt_[v][ip][1];
      pxdatum[v] = datum[v];
      trifocal3pt::invert_intrinsics(data::K_, pxdatum[v].col(ip).data(), datum[v].col(ip).data()); 
      trifocal3pt::invert_intrinsics_tgt(data::K_, pxdatum[v].col(ip).data()+2, datum[v].col(ip).data()+2);
    }
    {
    // 4th point is just the 2nd one, perturbed
    constexpr unsigned ip = 3;
    pxdatum[v].col(ip) = pxdatum[v].col(1);
    pxdatum[v](0,ip) += 5.0;
    pxdatum[v](1,ip) += 5.0;
    trifocal3pt::invert_intrinsics(data::K_, pxdatum[v].col(ip).data(), datum[v].col(ip).data()); 
    trifocal3pt::invert_intrinsics_tgt(data::K_, pxdatum[v].col(ip).data()+2, datum[v].col(ip).data()+2);
    }

    {
    // 5th and 6th are 1st repeated
    pxdatum[v].col(4) = pxdatum[v].col(0);
    datum[v].col(4) = datum[v].col(0);
    pxdatum[v].col(5) = pxdatum[v].col(0);
    datum[v].col(5) = datum[v].col(0);
    }
  }
  
  using TrifocalKernel = ThreeViewKernel<Trifocal3PointPositionTangentialSolver, 
                         Trifocal3PointPositionTangentialSolver>;
  
  const TrifocalKernel trifocal_kernel(datum[0], datum[1], datum[2], pxdatum[0], pxdatum[1], pxdatum[2], data::K_);
  
  double threshold = threshold_pixel_to_normalized(1.0, data::K_);
  threshold *= threshold; // squared error
  std::cerr << "Threshold squared in normalized coords "  << threshold << std::endl;
  unsigned constexpr max_iteration = 2; // testing
  // Vector of inliers for the best fit found
  vector<uint32_t> vec_inliers;
  const auto model = MaxConsensus(trifocal_kernel, 
      ScorerEvaluator<TrifocalKernel>(threshold), &vec_inliers, max_iteration);
  std::cerr << "Number of inliers (expect 5): "  << vec_inliers.size() << "\n";
  CHECK(vec_inliers.size() == 5);
  }
  
  { // 3 perfect points and n outliers
  }

  CHECK(true);
}

#if 0 // this test is ready, just needs some debugging
TEST(TrifocalSampleApp, fullrun) 
{
  TrifocalSampleApp T;
  int argcTest = 9;
  char argvTestData[9][50] = {"fullrun","-a" ,"/home/bielhpp/openMVG_bin/frame_00001.png","-b" ,"/home/bielhpp/openMVG_bin/frame_00030.png","-c" ,"/home/bielhpp/openMVG_bin/frame_00066.png","-K" ,"this is a stun"} ;
  char *argvTest[9];
  for (unsigned i = 0; i < argcTest ; i++) {
    argvTest[i] = (char *)argvTestData[i];
  }
  T.ProcessCmdLine(argcTest, argvTest);
  T.ExtractKeypoints();
  T.MatchKeypoints();
  T.ComputeTracks();
  //int n_ids = 5;
  //int desired_ids[n_ids] = {13, 23, 33, 63, 53};
  //CHECK(T.FilterIds(desired_ids, n_ids));
  T.Stats();
  T.ExtractXYOrientation();
  T.Display();
  T.DisplayDesiredIds();
  T.DisplayNonDesiredIds();
  T.RobustSolve();
  // T.DisplayInliers();
  T.DisplayInliersCamerasAndPoints();
  T.DisplayInliersCamerasAndPointsSIFT();
  std::cout<<"hej"<<std::endl;
  //return EXIT_SUCCESS;
  CHECK(true);
}
#endif

int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
