// Copyright (c) 2019 Pierre MOULON.
//:\file
//\Trifocal test with hardcoded samples
//\author Pierre MOULON
//\author Gabriel ANDRADE Rio de Janeiro State U.
//\author Ricardo Fabbri, Brown & Rio de Janeiro State U. (rfabbri.github.io) 

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <array>
#include <testing/testing.h>
#include "openMVG/multiview/trifocal/solver_trifocal_three_point.hpp"
#include "openMVG/multiview/trifocal/solver_trifocal_util.hpp"
#include "openMVG/cameras/Camera_Pinhole.hpp"

#include <minus/chicago-default.h>
#include <minus/internal-util.h>
#include <minus/debug-common.h>

using namespace std;
using namespace openMVG;
using namespace openMVG::trifocal;
using namespace openMVG::cameras;
typedef MiNuS::minus_util<double> util;

trifocal_model_t tt_gt_; // corresp. to minus' cameras_gt_

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
// - directly pass to solver
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
      Pinhole_Intrinsic::invert_intrinsics(data::K_, datum[v].col(ip).data(), datum[v].col(ip).data()); 
      Pinhole_Intrinsic::invert_intrinsics_tgt(data::K_, datum[v].col(ip).data()+2, datum[v].col(ip).data()+2);
    }
  }
  
  initialize_gt();
  // unsigned constexpr max_solve_tries = 1;
  bool found = false;

  // for (unsigned i = 0; i < max_solve_tries; ++i) {
    std::vector<trifocal_model_t> sols; // std::vector of 3-cam solutions
    
    // std::cerr << "Test log: Trying to solve, attempt: " << i+1 << std::endl;
    Trifocal3PointPositionTangentialSolver::Solve(datum[0], datum[1], datum[2], &sols);

    unsigned sol_id = (unsigned)-1;
    found = probe_solutions(sols, tt_gt_, &sol_id);
    if (found) {
      std::cerr << "Found solution at id " << sol_id << std::endl;
      //      for(unsigned j = 0; j < 3; j++)
      //        std::cout << sols[sol_id][j] << "\n" << std::endl;
      // break;
    }
    // std::cerr << "Test log: Solve failed to find ground truth. Retrying different randomization\n";
  // }
  
  CHECK(found);
}

int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
