// Copyright (c) 2019 Pierre MOULON.
//:\file
//\Trifocal test with hardcoded samples
//\author Pierre MOULON
//\author Gabriel ANDRADE Rio de Janeiro State U.
//\author Ricardo Fabbri, Brown & Rio de Janeiro State U. (rfabbri.github.io) 

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "trifocal-app.h"
#include "testing/testing.h"

using namespace trifocal3pt;
//-------------------------------------------------------------------------------
//
// Global variables

//---- HARDODED CASE SOLUTION -----------------------------------------------------
// Uncomment to hardcode solution for known case to test.
// This is as if MINUS was run and returned the following
// 
//constexpr unsigned n_ids = 5;
//unsigned desired_ids[n_ids] = {13, 23, 33, 93, 53};
//---- HARDODED CASE SOLUTION -----------------------------------------------------
// This is for hard ding test 
//class TrifocalAppTest : public TrifocalSampleApp {
//
//};
//
//TEST(TrifocalApp; Solve) {
//
//	double R0[3][3] = {
//		     {1,0,0},
//		     {0,1,0},
//		     {0,0,1}
//		    };
//
//	double T0[3][1] = {
//		       {13.022176},
//		       {-1.6546488},
//		       {352.47945}
//		      };
//
//	double R1[3][3] = {   
//		      {9.4083600000000001e-01,   2.7502399999999999e-01,   1.9796300000000000e-01},
//		      {2.9655799999999999e-01,  -3.8560499999999998e-01,  -8.7370599999999998e-01},
//		      {-1.6395499999999999e-01,   8.8072200000000000e-01,  -4.4435100000000000e-01}
//		     };
//	double T1[3][1] = { 
//		       {8.70420714},
//		       {-1.62157456},
//		       {-352.61248141}
//		     };
//
//	double R2[3][3] = {
//		      {0.970091125581631,   0.235130101826381,   0.060307903987350},
//		      {0.151694164781553,  -0.393265050435905,  -0.906824944780907},
//		      {-0.189504850701909,   0.888851188512892,  -0.417170799840706}
//		     
//	};
//	double T2[3][1] = { 
//		       {0.88920328},
//		       {-14.05063273},
//		       {-352.44248798}
//		     };
//
//	tt[0][0] = Mat34::Identity();
//	tt[0][1] = Mat34::Identity();
//	tt[0][2] = Mat34::Identity();
//	for(unsigned i=0;i<3;i++){
//		for(unsigned j=0;j<4;j++){
//			if(j<3){
//				tt[0][0](i,j) = R0[i][j];
//				tt[0][1](i,j) = R1[i][j];
//				tt[0][2](i,j) = R2[i][j];
//			}
//			else{
//				tt[0][0](i,j) = T0[i][1];
//				tt[0][1](i,j) = T1[i][1];
//				tt[0][2](i,j) = T2[i][1];
//			}
//		}                       
//	}
//	const solverResult=
//	const double delta = 0.001;     	
//	for(unsigned i=0;i<3;i++){
//		for(unsigned j=0;j<4;j++){
//			if(j<3){
//				EXPECT_NEAR(solverResult[0][0](i,j) , R0[i][j], );
//				EXPECT_NEAR(solverResult[0][1](i,j) , R1[i][j]);
//				EXPECT_NEAR(solverResult[0][2](i,j) , R2[i][j]);
//			}
//			else{
//				NEAR(solverResult[0][0](i,j) , T0[i][1]);
//				EXPECT_NEAR(solverResult[0][1](i,j) , T1[i][1]);
//				EXPECT_NEAR(solverResult[0][2](i,j) , T2[i][1]);
//			}
//		}                       
//	}
//}
//std::cerr << "Solutions:\n";
//for(unsigned i = 0;i < nsols_final; ++i){
//  for(unsigned j = 0;j < io::pp::nviews; ++j){
//    std::cerr << "Matrix" << i << " " << j << std::endl;
//    std::cerr << tt[i][j] << std::endl;
//    std::cerr << "\n";
//  }
//}
//  cout << "this is [R0|T0] " << "\n"; cout << tt[0][0] << "\n";
//  cout << "this is [R1|T1] " << "\n"; cout << tt[0][1] << "\n";
//  cout << "this is [R2|T2] " << "\n"; cout << tt[0][2] << "\n";


// Directly runs the solver and test
// - define synthetic data
// - directly passo to solver
//    - use Trifocal3PointPositionTangentialSolver::Solve
// - check if the solver returns any known root
// - might fail 5% of the time
// - also check if Trifocal3PointPositionTangentialSolver::error function returns zero
TEST(TrifocalSampleApp, solver) 
{
  trifocal_model_t tt; // std::vector of 3-cam solutions

  Mat d0, d1, d2;
  d0.resize(4,3);
  d1.resize(4,3);
  d2.resize(4,3);

  for (unsigned ip=0; ip < 3; ++ip) {
    d0(0,ip) = 
    d0(1,ip) = 
    d0(2,ip) = 
    d0(3,ip) = 
    
    d1(0,ip) = 
    d1(1,ip) = 
    d1(2,ip) = 
    d1(3,ip) = 
    
    d2(0,ip) = 
    d2(1,ip) = 
    d2(2,ip) = 
    d2(3,ip) = 
  }
  
  Trifocal3PointPositionTangentialSolver::Solve(d0, d1, d2, &tt);

  CHECK(true);
}

// Runs the solve through ransac 
// - first, synthetic data with three perfect points
// - second, synthetic data with one outlier
TEST(TrifocalSampleApp, solver) 
{
  CHECK(true);
}

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
  // T.DisplayInliersCamerasAndPointsSIFT();
  std::cout<<"hej"<<std::endl;
  //return EXIT_SUCCESS;
  CHECK(true);
}

int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
