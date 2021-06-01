
// Copyright (c) 2019 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/image/image_io.hpp"
#include "openMVG/image/image_concat.hpp"
#include "openMVG/features/akaze/image_describer_akaze.hpp"
#include "openMVG/features/sift/SIFT_Anatomy_Image_Describer.hpp"
#include "openMVG/matching/regions_matcher.hpp"
#include "openMVG/tracks/tracks.hpp"

#include "openMVG/numeric/extract_columns.hpp"
#include "openMVG/robust_estimation/robust_estimator_MaxConsensus.hpp"
//#include "openMVG/robust_estimation/robust_estimator_Ransac.hpp"
#include "openMVG/robust_estimation/score_evaluator.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "third_party/vectorGraphics/svgDrawer.hpp"

#include <memory>
#include <string>

//these are temporary includes, may be removed
#include <Eigen/StdVector>
#include <numeric>
#include "openMVG/multiview/projection.hpp"
#include "openMVG/multiview/triangulation.hpp"

#include "minus/minus.hxx"
#include "minus/chicago-default.h"
#include "minus/chicago14a-default-data.hxx"

#include "trifocal.h"

// Mat is Eigen::MatrixXd - matrix of doubles with dynamic size
// Vec3 is Eigen::Vector3d - Matrix< double, 3, 1 >

using namespace std;
using namespace openMVG;
using namespace openMVG::image;
using namespace openMVG::robust;
using namespace MiNuS;
using SIFT_Regions = openMVG::features::SIFT_Regions;


//Defining global variables

//constexpr unsigned n_ids = 5;
//unsigned desired_ids[n_ids] = {13, 23, 33, 93, 53};

//------------------------------------------------------------------------------
static void
revert_intrinsics(
    const double K[/*3 or 2 ignoring last line*/][3], 
    double pix_coords[2], 
    const double normalized_coords[2])
{
  double px = pix_coords;
  const double nrm = normalized_coords;
  // XXX: usar a inversa da formula exatamente como em invert_intrinsics.
  //      ter certeza que funciona se a entrada e saida forem mesmas posicoes de
  //      memoria
  px[0] = nrm[0]*K[0][0]+nrm[1]*K[0][1]+nrm[2]*K[0][2];
  px[1] = nrm[0]*K[1][0]+nrm[1]*K[1][1]+nrm[2]*K[1][2];
}

static void
revert_intrinsics_tgt(
    const double K[/*3 or 2 ignoring last line*/][3], 
    double pix_tgt_coords[2], 
    const double normalized_tgt_coords[2])
{
  double tp = pix_tgt_coords;
  const double t = normalized_tgt_coords;
  tp[0] = t[0]*K[0][0]+t[1]*K[0][1]+t[2]*K[0][2];
  tp[1] = t[0]*K[1][0]+t[1]*K[1][1]+t[2]*K[1][2];
}

static void
invert_intrinsics(
    const double K[/*3 or 2 ignoring last line*/][3], 
    const double pix_coords[2], 
    double normalized_coords[2])
{
  const double *px = pix_coords;
  double *nrm = normalized_coords;
  nrm[1] = (px[1] - K[1][2]) /K[1][1];
  nrm[0] = (px[0] - K[0][1]*nrm[1] - K[0][2])/K[0][0];
}

static void
invert_intrinsics_tgt(
    const double K[/*3 or 2 ignoring last line*/][3], 
    const double pix_tgt_coords[2], 
    double normalized_tgt_coords[2])
{
  const double *tp = pix_tgt_coords;
  double *t = normalized_tgt_coords;
  t[1] = tp[1]/K[1][1];
  t[0] = (tp[0] - K[0][1]*tp[1])/K[0][0];
}
// See big notes eq. 5.2.13 at beginning of the code.


//-------------------------------------------------------------------------------
//
// Global variables
// 
// constexpr unsigned n_ids = 5;
// unsigned desired_ids[n_ids] = {13, 23, 33, 93, 53};
//-------------------------------------------------------------------------------


//-------------------------------------------------------------------------------
// Trifocal3PointPositionTangentialSolver
//-------------------------------------------------------------------------------

void Trifocal3PointPositionTangentialSolver::
Solve(
      const Mat &datum_0,
      const Mat &datum_1,
      const Mat &datum_2,
      std::vector<trifocal_model_t> *trifocal_tensor)
{
  double p[io::pp::nviews][io::pp::npoints][io::ncoords2d];
  double tgt[io::pp::nviews][io::pp::npoints][io::ncoords2d]; 

  std::cerr << "TRIFOCAL LOG: Called Solve()\n";
  // pack into solver's efficient representation
  for (unsigned ip=0; ip < io::pp::npoints; ++ip) {
    p[0][ip][0] = datum_0(0,ip);
    p[0][ip][1] = datum_0(1,ip);
    tgt[0][ip][0] = datum_0(2,ip); 
    tgt[0][ip][1] = datum_0(3,ip); 
    
    p[1][ip][0] = datum_1(0,ip);
    p[1][ip][1] = datum_1(1,ip);
    tgt[1][ip][0] = datum_1(2,ip); 
    tgt[1][ip][1] = datum_1(3,ip); 
    
    p[2][ip][0] = datum_2(0,ip);
    p[2][ip][1] = datum_2(1,ip);
    tgt[2][ip][0] = datum_2(2,ip); 
    tgt[2][ip][1] = datum_2(3,ip); 
  }
  
  unsigned nsols_final = 0;
  unsigned id_sols[M::nsols];
  double  cameras[M::nsols][io::pp::nviews-1][4][3];  // first camera is always [I | 0]
  
  std::cerr << "TRIFOCAL LOG: Before minus::solve()\n" << std::endl;
  MiNuS::minus<chicago>::solve(p, tgt, cameras, id_sols, &nsols_final);
  //std::cerr << datum_0  "\n"; 
    
  // double R0[3][3] = {
  //                    {1,0,0},
  //                    {0,1,0},
  //                    {0,0,1}
  //                   };
  //   
  //  double T0[3][1] = {
  //                      {13.022176},
  //                      {-1.6546488},
  //                      {352.47945}
  //                     };

  //  double R1[3][3] = {   
  //                     {9.4083600000000001e-01,   2.7502399999999999e-01,   1.9796300000000000e-01},
  //                     {2.9655799999999999e-01,  -3.8560499999999998e-01,  -8.7370599999999998e-01},
  //                     {-1.6395499999999999e-01,   8.8072200000000000e-01,  -4.4435100000000000e-01}
  //                    };
  //  double T1[3][1] = { 
  //                      {8.70420714},
  //                      {-1.62157456},
  //                      {-352.61248141}
  //                    };

  //  double R2[3][3] = {
  //                     {0.970091125581631,   0.235130101826381,   0.060307903987350},
  //                     {0.151694164781553,  -0.393265050435905,  -0.906824944780907},
  //                     {-0.189504850701909,   0.888851188512892,  -0.417170799840706}
  //                    };

  //  double T2[3][1] = { 
  //                      {0.88920328},
  //                      {-14.05063273},
  //                      {-352.44248798}
  //                    };
  
 //fill C0* with for loop
  std::cerr << "Number of sols " << nsols_final << std::endl;
  std::vector<trifocal_model_t> &tt = *trifocal_tensor; // if I use the STL container, 
  // This I would have to change the some other pieces of code, maybe altering the entire logic of this program!!
  // std::cerr << "TRIFOCAL LOG: Antes de resize()\n" << std::endl;
  tt.resize(nsols_final);
  std::cerr << "TRIFOCAL LOG: Chamou resize()\n";
  //using trifocal_model_t = array<Mat34, 3>;
  for (unsigned s=0; s < nsols_final; ++s) {
    tt[s][0] = Mat34::Identity(); // view 0 [I | 0]
    for (unsigned v=1; v < io::pp::nviews; ++v) {
        memcpy(tt[s][v].data(), (double *) cameras[id_sols[s]][v], 9*sizeof(double));
        for (unsigned r=0; r < 3; ++r)
          tt[s][v](r,3) = cameras[id_sols[s]][v][3][r];
    }
  }
  //This is for hard coding test 
  //tt[0][0] = Mat34::Identity();
  //tt[0][1] = Mat34::Identity();
  //tt[0][2] = Mat34::Identity();
  //for(unsigned i=0;i<3;i++){
  //  for(unsigned j=0;j<4;j++){
  //    if(j<3){
  //      tt[0][0](i,j) = R0[i][j];
  //      tt[0][1](i,j) = R1[i][j];
  //      tt[0][2](i,j) = R2[i][j];
  //    }
  //    else{
  //      tt[0][0](i,j) = T0[i][1];
  //      tt[0][1](i,j) = T1[i][1];
  //      tt[0][2](i,j) = T2[i][1];
  //    }
  //  }                       
  //}
  // cout << "this is [R0|T0] " << "\n"; cout << tt[0][0] << "\n";
  // cout << "this is [R1|T1] " << "\n"; cout << tt[0][1] << "\n";
  // cout << "this is [R2|T2] " << "\n"; cout << tt[0][2] << "\n";
  
  // TODO: filter the solutions by:
  // - positive depth and 
  // - using tangent at 3rd point
  //
  //  if we know the rays are perfectly coplanar, we can just use cross
  // product within the plane instead of SVD
  std::cerr << "TRIFOCAL LOG: Finished ()Solve()\n";
}

void Trifocal3PointPositionTangentialSolver::
static double Error(
  const trifocal_model_t &tt,
  const Vec &bearing_0, // x,y,tangentialx,tangentialy
  const Vec &bearing_1,
  const Vec &bearing_2,
  const Vec &pixbearing_0,
  const Vec &pixbearing_1,
  const Vec &pixbearing_2,
  const double K_[2][3]) 
{
  //std::cerr << "TRIFOCAL LOG: Called Error()\n";
  // Return the cost related to this model and those sample data point
  // Ideal algorithm:
  // 1) reconstruct the 3D points and orientations
  // 2) project the 3D points and orientations on all images_
  // 3) compute error 
  // 
  // In practice we ignore the directions and only reproject to one third view
  
  // 3x3: each column is x,y,1
  Mat3 bearing;
  bearing << bearing_0.head(2).homogeneous(),
             bearing_1.head(2).homogeneous(), 
             bearing_2.head(2).homogeneous();
  Mat2 pixbearing; // << XXX mat2
  pixbearing << pixbearing_0.head(2).homogeneous(),
                pixbearing_1.head(2).homogeneous();
  // Using triangulation.hpp
  Vec4 triangulated_homg;
  unsigned third_view = 0;
  // pick the wider baseline. TODO: measure all pairwise translation distances
  if (tt[1].col(3).squaredNorm() > tt[2].col(3).squaredNorm()) {
    // TODO use triangulation from the three views at once
    TriangulateDLT(tt[0], bearing.col(0), tt[1], bearing.col(1), &triangulated_homg);
    third_view = 2;
  } else {
    TriangulateDLT(tt[0], bearing.col(0), tt[2], bearing.col(2), &triangulated_homg);
    third_view = 1;
  }
  
  // Computing the projection of triangulated points using projection.hpp
  // For prototyping and speed, for now we will only project to the third view
  // and report only one error
  Vec2 pxreprojected = Vec3(tt[third_view]*triangulated_homg).hnormalized();
  // XXX revert intrinsics to measure the error in pixels
  revert_intrinsics(K, pxreprojected, pxreprojected);
   
  Vec2 pxmeasured    = pxbearing.col(third_view);
  //cout << "error " << (reprojected - measured).squaredNorm() << "\n";
  //cout << "triang " <<triangulated_homg <<"\n";
  //std::cerr << "TRIFOCAL LOG: Finished Error()\n";
  return (reprojected-measured).squaredNorm();
}

//-------------------------------------------------------------------------------
static int iteration_global_debug = 0;


int main(int argc, char **argv) {
  TrifocalSampleApp T;
  
  T.ProcessCmdLine(argc, argv);
  T.ExtractKeypoints();
  T.MatchKeypoints();
  T.ComputeTracks();
  T.Stats();
  T.ExtractXYOrientation();
  T.Display();
  T.DisplayDesiredIds();
  T.RobustSolve();
  // T.DisplayInliers();
  T.DisplayInliersCamerasAndPoints();
  // T.DisplayInliersCamerasAndPointsSIFT();

  return EXIT_SUCCESS;
}
