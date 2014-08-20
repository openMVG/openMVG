
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include "openMVG/cameras/PinholeCamera.hpp"
#include "openMVG/image/image.hpp"
#include "openMVG/features/features.hpp"
#include "openMVG/matching/matcher_brute_force.hpp"
#include "openMVG/matching/indMatchDecoratorXY.hpp"
#include "openMVG/multiview/projection.hpp"
#include "openMVG/multiview/triangulation.hpp"
#include "openMVG_Samples/robust_essential/essential_estimation.hpp"

#include "openMVG_Samples/siftPutativeMatches/two_view_matches.hpp"

// Bundle Adjustment includes
#include "openMVG/bundle_adjustment/problem_data_container.hpp"
#include "openMVG/bundle_adjustment/pinhole_ceres_functor.hpp"
#include "openMVG/bundle_adjustment/pinhole_brown_Rt_ceres_functor.hpp"

#include "patented/sift/SIFT.hpp"

#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "third_party/vectorGraphics/svgDrawer.hpp"

#include <string>
#include <iostream>

using namespace openMVG;
using namespace openMVG::matching;
using namespace openMVG::bundle_adjustment;
using namespace svg;
using namespace std;

/// Read intrinsic K matrix from a file (ASCII)
/// F 0 ppx
/// 0 F ppy
/// 0 0 1
bool readIntrinsic(const std::string & fileName, Mat3 & K);

/// Export 3D point vector and camera position to PLY format
bool exportToPly(const std::vector<Vec3> & vec_points,
  const std::vector<Vec3> & vec_camPos,
  const std::string & sFileName);

/// Triangulate and export valid point as PLY (point in front of the cameras)
void triangulateAndSaveResult(
  const PinholeCamera & camL,
  const PinholeCamera & camR,
  std::vector<size_t> & vec_inliers,
  const Mat & xL,
  const Mat & xR,
  std::vector<Vec3> & vec_3DPoints);

/// Perform a Bundle Adjustment: Refine the camera [R|t|focal] and the structure
void do_bundle_adjustment(
  PinholeCamera & camL,
  PinholeCamera & camR,
  const Mat & xL,
  const Mat & xR,
  const std::vector<size_t> & vec_inliers,
  std::vector<Vec3> & vec_3DPoints);

/// Perform a Bundle Adjustment: Refine the cameras [R|t],
///  and common pinhole intrinsic [focal,ppx,ppy] and the structure
void do_bundle_adjustment_common_intrinsic_pinhole(
  PinholeCamera & camL,
  PinholeCamera & camR,
  const Mat & xL,
  const Mat & xR,
  const std::vector<size_t> & vec_inliers,
  std::vector<Vec3> & vec_3DPoints);

/// Perform a Bundle Adjustment: Refine the cameras [R|t]
///  and common intrinsics [focal,ppx,ppy,k1,k2,k3] and the structure
void do_bundle_adjustment_common_intrinsics_brown_pinhole(
  PinholeCamera & camL,
  PinholeCamera & camR,
  const Mat & xL,
  const Mat & xR,
  const std::vector<size_t> & vec_inliers,
  std::vector<Vec3> & vec_3DPoints);

/// Show :
///  how computing an essential with know internal calibration matrix K
///  how refine the camera motion, focal and structure with Bundle Adjustment
///   way 1: independent cameras [R|t|f] and structure
///   way 2: independent cameras motion [R|t], shared focal [f] and structure
int main() {

  std::string sInputDir = stlplus::folder_up(string(THIS_SOURCE_DIR))
    + "/imageData/SceauxCastle/";
  Image<RGBColor> image;
  string jpg_filenameL = sInputDir + "100_7101.jpg";
  string jpg_filenameR = sInputDir + "100_7102.jpg";

  Image<unsigned char> imageL, imageR;
  ReadImage(jpg_filenameL.c_str(), &imageL);
  ReadImage(jpg_filenameR.c_str(), &imageR);

  // Define the used descriptor (SIFT : 128 float value)
  typedef float descType;
  typedef Descriptor<descType,128> SIFTDescriptor;

  // Prepare vector to store detected feature and associated descriptor
  vector<SIOPointFeature> featsL, featsR;
  vector<SIFTDescriptor > descsL, descsR;
  // Call SIFT detector
  bool bOctaveMinus1 = false;
  bool bRootSift = true;
  SIFTDetector(imageL, featsL, descsL, bOctaveMinus1, bRootSift);
  SIFTDetector(imageR, featsR, descsR, bOctaveMinus1, bRootSift);

  // Show both images side by side
  {
    Image<unsigned char> concat;
    ConcatH(imageL, imageR, concat);
    string out_filename = "01_concat.jpg";
    WriteImage(out_filename.c_str(), concat);
  }

  //- Draw features on the two image (side by side)
  {
    Image<unsigned char> concat;
    ConcatH(imageL, imageR, concat);

    //-- Draw features :
    for (size_t i=0; i < featsL.size(); ++i )  {
      const SIOPointFeature & imaA = featsL[i];
      DrawCircle(imaA.x(), imaA.y(), imaA.scale(), 255, &concat);
    }
    for (size_t i=0; i < featsR.size(); ++i )  {
      const SIOPointFeature & imaB = featsR[i];
      DrawCircle(imaB.x()+imageL.Width(), imaB.y(), imaB.scale(), 255, &concat);
    }
    string out_filename = "02_features.jpg";
    WriteImage(out_filename.c_str(), concat);
  }

  std::vector<IndMatch> vec_PutativeMatches;
  //-- Perform matching -> find Nearest neighbor, filtered with Distance ratio
  {
    // Define the matcher
    //  and the used metric (Squared L2)
    typedef L2_Vectorized<SIFTDescriptor::bin_type> Metric;
    // Brute force matcher is defined as following:
    typedef ArrayMatcherBruteForce<SIFTDescriptor::bin_type, Metric> MatcherT;

    // Distance ratio quite high in order to have noise corrupted data. Squared due to squared metric
    getPutativesMatches<SIFTDescriptor, MatcherT>(descsL, descsR, Square(0.8), vec_PutativeMatches);

    IndMatchDecorator<float> matchDeduplicator(
            vec_PutativeMatches, featsL, featsR);
    matchDeduplicator.getDeduplicated(vec_PutativeMatches);

    // Draw correspondences after Nearest Neighbor ratio filter
    svgDrawer svgStream( imageL.Width() + imageR.Width(), max(imageL.Height(), imageR.Height()));
    svgStream.drawImage(jpg_filenameL, imageL.Width(), imageL.Height());
    svgStream.drawImage(jpg_filenameR, imageR.Width(), imageR.Height(), imageL.Width());
    for (size_t i = 0; i < vec_PutativeMatches.size(); ++i) {
      //Get back linked feature, draw a circle and link them by a line
      const SIOPointFeature & L = featsL[vec_PutativeMatches[i]._i];
      const SIOPointFeature & R = featsR[vec_PutativeMatches[i]._j];
      svgStream.drawLine(L.x(), L.y(), R.x()+imageL.Width(), R.y(), svgStyle().stroke("green", 2.0));
      svgStream.drawCircle(L.x(), L.y(), L.scale(), svgStyle().stroke("yellow", 2.0));
      svgStream.drawCircle(R.x()+imageL.Width(), R.y(), R.scale(),svgStyle().stroke("yellow", 2.0));
    }
    string out_filename = "03_siftMatches.svg";
    ofstream svgFile( out_filename.c_str() );
    svgFile << svgStream.closeSvgFile().str();
    svgFile.close();
  }

  // Essential geometry filtering of putative matches
  {
    Mat3 K;
    //read K from file
    if (!readIntrinsic(stlplus::create_filespec(sInputDir,"K","txt"), K))
    {
      std::cerr << "Cannot read intrinsic parameters." << std::endl;
      return EXIT_FAILURE;
    }

    //A. prepare the corresponding putatives points
    Mat xL(2, vec_PutativeMatches.size());
    Mat xR(2, vec_PutativeMatches.size());
    for (size_t k = 0; k < vec_PutativeMatches.size(); ++k)  {
      const SIOPointFeature & imaL = featsL[vec_PutativeMatches[k]._i];
      const SIOPointFeature & imaR = featsR[vec_PutativeMatches[k]._j];
      xL.col(k) = imaL.coords().cast<double>();
      xR.col(k) = imaR.coords().cast<double>();
    }

    //B. robust estimation of the essential matrix
    std::vector<size_t> vec_inliers;
    Mat3 E;
    std::pair<size_t, size_t> size_imaL(imageL.Width(), imageL.Height());
    std::pair<size_t, size_t> size_imaR(imageR.Width(), imageR.Height());
    double thresholdE = 0.0, NFA = 0.0;
    if (robustEssential(
      K, K,         // intrinsics
      xL, xR,       // corresponding points
      &E,           // essential matrix
      &vec_inliers, // inliers
      size_imaL,    // Left image size
      size_imaR,    // Right image size
      &thresholdE,  // Found AContrario Theshold
      &NFA,         // Found AContrario NFA
      std::numeric_limits<double>::infinity()))
    {
      std::cout << "\nFound an Essential matrix under the confidence threshold of: "
        << thresholdE << " pixels\n\twith: " << vec_inliers.size() << " inliers"
        << " from: " << vec_PutativeMatches.size()
        << " putatives correspondences"
        << std::endl;

      // Show Essential validated point
      svgDrawer svgStream( imageL.Width() + imageR.Width(), max(imageL.Height(), imageR.Height()));
      svgStream.drawImage(jpg_filenameL, imageL.Width(), imageL.Height());
      svgStream.drawImage(jpg_filenameR, imageR.Width(), imageR.Height(), imageL.Width());
      for ( size_t i = 0; i < vec_inliers.size(); ++i)  {
        const SIOPointFeature & LL = featsL[vec_PutativeMatches[vec_inliers[i]]._i];
        const SIOPointFeature & RR = featsR[vec_PutativeMatches[vec_inliers[i]]._j];
        const Vec2f L = LL.coords();
        const Vec2f R = RR.coords();
        svgStream.drawLine(L.x(), L.y(), R.x()+imageL.Width(), R.y(), svgStyle().stroke("green", 2.0));
        svgStream.drawCircle(L.x(), L.y(), LL.scale(), svgStyle().stroke("yellow", 2.0));
        svgStream.drawCircle(R.x()+imageL.Width(), R.y(), RR.scale(),svgStyle().stroke("yellow", 2.0));
      }
      string out_filename = "04_ACRansacEssential.svg";
      ofstream svgFile( out_filename.c_str() );
      svgFile << svgStream.closeSvgFile().str();
      svgFile.close();

      //C. Extract the rotation and translation of the camera from the essential matrix
      Mat3 R;
      Vec3 t;
      if (!estimate_Rt_fromE(K, K, xL, xR, E, vec_inliers,
        &R, &t))
      {
        std::cerr << " /!\\ Failed to compute initial R|t for the initial pair" << std::endl;
        return false;
      }
      std::cout << std::endl
        << "-- Rotation|Translation matrices: --" << std::endl
        << R << std::endl << std::endl << t << std::endl;

      // Build Left and Right Camera
      PinholeCamera camL(K, Mat3::Identity(), Vec3::Zero());
      PinholeCamera camR(K, R, t);

      //C. Triangulate and check valid points
      // invalid point that do not respect cheirality are discarded (removed
      //  from the list of inliers.
      std::vector<Vec3> vec_3DPoints;
      triangulateAndSaveResult(
        camL, camR,
        vec_inliers,
        xL, xR, vec_3DPoints);

      //D. Refine the computed structure and cameras
      std::cout << "Which BA do you want ? " << std::endl
        << "\t 1: Refine [X],[f,R|t] (individual cameras)\n"
        << "\t 2: Refine [X],[R|t], shared [f, ppx, ppy]\n"
        << "\t 3: Refine [X],[R|t], shared brown disto models [f,ppx,ppy,k1,k2,k3]\n" << std::endl;
      int iBAType = -1;
      std::cin >> iBAType;
      switch(iBAType)
      {
        case 1:
        {
          do_bundle_adjustment(
            camL, camR,
            xL, xR,
            vec_inliers,
            vec_3DPoints);
        }
        break;

        case 2:
        {
          do_bundle_adjustment_common_intrinsic_pinhole(
            camL, camR,
            xL, xR,
            vec_inliers,
            vec_3DPoints);
        }
        break;

        case 3:
        {
          do_bundle_adjustment_common_intrinsics_brown_pinhole(
            camL, camR,
            xL, xR,
            vec_inliers,
            vec_3DPoints);
        }
        break;
        default:
          std::cerr << "Invalid input number" << std::endl;
      }


      //E. Export as PLY the refined scene (camera pos + 3Dpoints)
      std::vector<Vec3> vec_camPos;
      vec_camPos.push_back( camL._C );
      vec_camPos.push_back( camR._C );
      exportToPly(vec_3DPoints, vec_camPos, "EssentialGeometry.ply");

    }
    else  {
      std::cout << "ACRANSAC was unable to estimate a rigid essential matrix"
        << std::endl;
    }
  }
  return EXIT_SUCCESS;
}

bool readIntrinsic(const std::string & fileName, Mat3 & K)
{
  // Load the K matrix
  ifstream in;
  in.open( fileName.c_str(), ifstream::in);
  if(in.is_open())  {
    for (int j=0; j < 3; ++j)
      for (int i=0; i < 3; ++i)
        in >> K(j,i);
  }
  else  {
    std::cerr << std::endl
      << "Invalid input K.txt file" << std::endl;
    return false;
  }
  return true;
}

/// Export 3D point vector and camera position to PLY format
bool exportToPly(const std::vector<Vec3> & vec_points,
  const std::vector<Vec3> & vec_camPos,
  const std::string & sFileName)
{
  std::ofstream outfile;
  outfile.open(sFileName.c_str(), std::ios_base::out);

  outfile << "ply"
    << '\n' << "format ascii 1.0"
    << '\n' << "element vertex " << vec_points.size()+vec_camPos.size()
    << '\n' << "property float x"
    << '\n' << "property float y"
    << '\n' << "property float z"
    << '\n' << "property uchar red"
    << '\n' << "property uchar green"
    << '\n' << "property uchar blue"
    << '\n' << "end_header" << std::endl;

  for (size_t i=0; i < vec_points.size(); ++i)  {
      outfile << vec_points[i].transpose()
      << " 255 255 255" << "\n";
  }

  for (size_t i=0; i < vec_camPos.size(); ++i)  {
    outfile << vec_camPos[i].transpose()
      << " 0 255 0" << "\n";
  }
  outfile.flush();
  bool bOk = outfile.good();
  outfile.close();
  return bOk;
}

/// Triangulate and export valid point as PLY (point in front of the cameras)
void triangulateAndSaveResult(
  const PinholeCamera & camL,
  const PinholeCamera & camR,
  std::vector<size_t> & vec_inliers,
  const Mat & xL,
  const Mat & xR,
  std::vector<Vec3> & vec_3DPoints)
{
  size_t nb_invalid3DPoints = 0;
  std::vector<size_t> vec_valid3DPoints;
  std::vector<double> vec_residuals;
  for (size_t k = 0; k < vec_inliers.size(); ++k) {
    const Vec2 & xL_ = xL.col(vec_inliers[k]);
    const Vec2 & xR_ = xR.col(vec_inliers[k]);

    Vec3 X = Vec3::Zero();
    TriangulateDLT(camL._P, xL_, camR._P, xR_, &X);

    // Compute residual:
    double dResidual = (camL.Residual(X, xL_) + camR.Residual(X, xR_))/2.0;
    vec_residuals.push_back(dResidual);
    if (camL.Depth(X) < 0 && camR.Depth(X) < 0) {
      ++nb_invalid3DPoints;
    }
    else  {
      vec_3DPoints.push_back(X);
      vec_valid3DPoints.push_back(vec_inliers[k]);
    }
  }
  if (nb_invalid3DPoints > 0)
  {
    std::cout << nb_invalid3DPoints
      << " correspondence(s) with negative depth have been discarded."
      << std::endl;
    // remove from the inliers list the point that are behind the camera
    vec_inliers = vec_valid3DPoints;
  }

  // Display some statistics of reprojection errors
  float dMin, dMax, dMean, dMedian;
  minMaxMeanMedian<float>(vec_residuals.begin(), vec_residuals.end(),
                        dMin, dMax, dMean, dMedian);

  std::cout << std::endl
    << "Essential matrix estimation, residuals statistics:" << "\n"
    << "\t-- Residual min:\t" << dMin << std::endl
    << "\t-- Residual median:\t" << dMedian << std::endl
    << "\t-- Residual max:\t "  << dMax << std::endl
    << "\t-- Residual mean:\t " << dMean << std::endl;
}

void do_bundle_adjustment(
  PinholeCamera & camL,
  PinholeCamera & camR,
  const Mat & xL,
  const Mat & xR,
  const std::vector<size_t> & vec_inliers,
  std::vector<Vec3> & vec_3DPoints)
{
  int nviews = 2;
  int n3Dpoints = vec_inliers.size();

  // Setup a BA problem
  BA_Problem_data<7> ba_problem;

  // Configure the size of the problem
  ba_problem.num_cameras_ = nviews;
  ba_problem.num_points_ = n3Dpoints;
  ba_problem.num_observations_ = nviews * n3Dpoints;

  ba_problem.point_index_.reserve(ba_problem.num_observations_);
  ba_problem.camera_index_.reserve(ba_problem.num_observations_);
  ba_problem.observations_.reserve(2 * ba_problem.num_observations_);

  ba_problem.num_parameters_ = 7 * ba_problem.num_cameras_ + 3 * ba_problem.num_points_;
  ba_problem.parameters_.reserve(ba_problem.num_parameters_);

  // Fill it with data (For each 3D point setup the tracks : the 2D visbility)
  PinholeCamera vec_cam[2] = {camL, camR};
  for (int i = 0; i < n3Dpoints; ++i) {
    // Collect the image of point i in each frame (xL, xR).
    const Vec2 & xL_ = xL.col(vec_inliers[i]);
    const Vec2 & xR_ = xR.col(vec_inliers[i]);

    // Left 2D observations
    double ppx = vec_cam[0]._K(0,2), ppy = vec_cam[0]._K(1,2);
    ba_problem.camera_index_.push_back(0);
    ba_problem.point_index_.push_back(i);
    ba_problem.observations_.push_back( xL_(0) - ppx);
    ba_problem.observations_.push_back( xL_(1) - ppy);

    // Right 2D observations
    ppx = vec_cam[1]._K(0,2);
    ppy = vec_cam[1]._K(1,2);
    ba_problem.camera_index_.push_back(1);
    ba_problem.point_index_.push_back(i);
    ba_problem.observations_.push_back( xR_(0) - ppx);
    ba_problem.observations_.push_back( xR_(1) - ppy);
  }

  // Add camera parameters (R, t, focal)
  for (int j = 0; j < nviews; ++j) {
    // Rotation matrix to angle axis
    std::vector<double> angleAxis(3);
    ceres::RotationMatrixToAngleAxis((const double*)vec_cam[j]._R.data(), &angleAxis[0]);
    // translation
    Vec3 t = vec_cam[j]._t;
    double focal = vec_cam[j]._K(0,0);
    ba_problem.parameters_.push_back(angleAxis[0]);
    ba_problem.parameters_.push_back(angleAxis[1]);
    ba_problem.parameters_.push_back(angleAxis[2]);
    ba_problem.parameters_.push_back(t[0]);
    ba_problem.parameters_.push_back(t[1]);
    ba_problem.parameters_.push_back(t[2]);
    ba_problem.parameters_.push_back(focal);
  }

  // Add 3D points coordinates parameters
  for (int i = 0; i < n3Dpoints; ++i) {
    Vec3 pt3D = vec_3DPoints[i];
    ba_problem.parameters_.push_back(pt3D[0]);
    ba_problem.parameters_.push_back(pt3D[1]);
    ba_problem.parameters_.push_back(pt3D[2]);
  }

  // Create residuals for each observation in the bundle adjustment problem. The
  // parameters for cameras and points are added automatically.
  ceres::Problem problem;
  for (int i = 0; i < ba_problem.num_observations(); ++i) {

    // Each Residual block takes a point and a camera as input and outputs a 2
    // dimensional residual. Internally, the cost function stores the observed
    // image location and compares the reprojection against the observation.
    ceres::CostFunction* cost_function =
        new ceres::AutoDiffCostFunction<pinhole_reprojectionError::ErrorFunc_Refine_Camera_3DPoints, 2, 7, 3>(
            new pinhole_reprojectionError::ErrorFunc_Refine_Camera_3DPoints(
                & ba_problem.observations()[2 * i]));

    problem.AddResidualBlock(cost_function,
                             NULL, // squared loss
                             ba_problem.mutable_camera_for_observation(i),
                             ba_problem.mutable_point_for_observation(i));
  }

  // Make Ceres automatically detect the bundle structure. Note that the
  // standard solver, SPARSE_NORMAL_CHOLESKY, also works fine but it is slower
  // for standard bundle adjustment problems.
  ceres::Solver::Options options;
  options.linear_solver_type = ceres::SPARSE_SCHUR;
  if (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::SUITE_SPARSE))
    options.sparse_linear_algebra_library_type = ceres::SUITE_SPARSE;
  else
    if (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::CX_SPARSE))
      options.sparse_linear_algebra_library_type = ceres::CX_SPARSE;
    else
    {
      // No sparse backend for Ceres.
      // Use dense solving
      options.linear_solver_type = ceres::DENSE_SCHUR;
    }
  options.minimizer_progress_to_stdout = false;
  options.logging_type = ceres::SILENT;

  ceres::Solver::Summary summary;
  ceres::Solve(options, &problem, &summary);

  double dResidual_before = std::sqrt( summary.initial_cost / (ba_problem.num_observations_*2.));
  double dResidual_after = std::sqrt( summary.final_cost / (ba_problem.num_observations_*2.));

  std::cout << std::endl
    << "Bundle Adjustment of cameras [R|t|f] and Structure : \n"
    << " Initial RMSE : " << dResidual_before << "\n"
    << " Final RMSE : " << dResidual_after << std::endl;

  // If no error, get back refined parameters
  if (summary.IsSolutionUsable())
  {
    // Get back 3D points
    size_t cpt = 0;
    for (std::vector<Vec3>::iterator iter = vec_3DPoints.begin();
      iter != vec_3DPoints.end(); ++iter, ++cpt)
    {
      const double * pt = ba_problem.mutable_points() + cpt*3;
      Vec3 & pt3D = *iter;
      pt3D = Vec3(pt[0], pt[1], pt[2]);
    }
    // Get back camera
    for (cpt = 0; cpt < nviews; ++cpt)
    {
      const double * cam = ba_problem.mutable_cameras() + cpt*7;
      Mat3 R;
      // angle axis to rotation matrix
      ceres::AngleAxisToRotationMatrix(cam, R.data());
      Vec3 t(cam[3], cam[4], cam[5]);
      double focal = cam[6];

      // Update the camera
      PinholeCamera & sCam = vec_cam[cpt];
      Mat3 K = sCam._K;
      K(0,0) = K(1,1) = focal;
      std::cout << "Refined focal[" << cpt << "]: " << focal << std::endl;
      sCam = PinholeCamera(K, R, t);
    }
  }
}

/// Perform a Bundle Adjustment: Refine the cameras [R|t],
///  common intrinsic [focal,ppx,ppy] and the structure
void do_bundle_adjustment_common_intrinsic_pinhole(
  PinholeCamera & camL,
  PinholeCamera & camR,
  const Mat & xL,
  const Mat & xR,
  const std::vector<size_t> & vec_inliers,
  std::vector<Vec3> & vec_3DPoints)
{
  int nCameraMotion = 2;
  int nCameraIntrinsic = 1;
  int n3Dpoints = vec_inliers.size();

  // Setup a BA problem
  BA_Problem_data_camMotionAndIntrinsic<6,3> ba_problem;

  // Configure the size of the problem
  ba_problem.num_cameras_ = nCameraMotion;
  ba_problem.num_intrinsics_ = nCameraIntrinsic;
  ba_problem.num_points_ = n3Dpoints;
  ba_problem.num_observations_ = nCameraMotion * n3Dpoints;

  ba_problem.point_index_.reserve(ba_problem.num_observations_);
  ba_problem.camera_index_extrinsic.reserve(ba_problem.num_observations_);
  ba_problem.camera_index_intrinsic.reserve(ba_problem.num_observations_);
  ba_problem.observations_.reserve(2 * ba_problem.num_observations_);

  ba_problem.num_parameters_ =
    6 * ba_problem.num_cameras_ + 3 * ba_problem.num_intrinsics_ + 3 * ba_problem.num_points_;
  ba_problem.parameters_.reserve(ba_problem.num_parameters_);

  // Fill it with data (For each 3D point setup the tracks : the 2D visbility)
  // The two camera share the same intrinsic
  PinholeCamera vec_cam[2] = {camL, camR};
  for (int i = 0; i < n3Dpoints; ++i) {
    // Collect the image of point i in each frame (xL, xR).
    const Vec2 & xL_ = xL.col(vec_inliers[i]);
    const Vec2 & xR_ = xR.col(vec_inliers[i]);

    // Left 2D observations

    ba_problem.camera_index_extrinsic.push_back(0);
    ba_problem.camera_index_intrinsic.push_back(0);
    ba_problem.point_index_.push_back(i);
    ba_problem.observations_.push_back(xL_(0));
    ba_problem.observations_.push_back(xL_(1));

    // Right 2D observations
    ba_problem.camera_index_extrinsic.push_back(1);
    ba_problem.camera_index_intrinsic.push_back(0); // same intrinsic group
    ba_problem.point_index_.push_back(i);
    ba_problem.observations_.push_back(xR_(0));
    ba_problem.observations_.push_back(xR_(1));
  }

  // Add camera extrinsics [R,t]
  for (int j = 0; j < nCameraMotion; ++j) {
    // Rotation matrix to angle axis
    std::vector<double> angleAxis(3);
    ceres::RotationMatrixToAngleAxis((const double*)vec_cam[j]._R.data(), &angleAxis[0]);
    // translation
    Vec3 t = vec_cam[j]._t;
    ba_problem.parameters_.push_back(angleAxis[0]);
    ba_problem.parameters_.push_back(angleAxis[1]);
    ba_problem.parameters_.push_back(angleAxis[2]);
    ba_problem.parameters_.push_back(t[0]);
    ba_problem.parameters_.push_back(t[1]);
    ba_problem.parameters_.push_back(t[2]);
  }
  // Add camera intrinsic (focal)
  double
    focal = (vec_cam[0]._K(0,0) + vec_cam[0]._K(1,1)
     + vec_cam[1]._K(1,1) + vec_cam[1]._K(0,0)) / 4.0;
  double ppx = (vec_cam[0]._K(0,2) + vec_cam[1]._K(0,2)) / 2.0;
  double ppy = (vec_cam[0]._K(1,2) + vec_cam[1]._K(1,2)) / 2.0;
  // Setup the intrinsic in the ba_problem
  ba_problem.parameters_.push_back(focal);
  ba_problem.parameters_.push_back(ppx);
  ba_problem.parameters_.push_back(ppy);

  // Add 3D points coordinates parameters
  for (int i = 0; i < n3Dpoints; ++i) {
    Vec3 pt3D = vec_3DPoints[i];
    ba_problem.parameters_.push_back(pt3D[0]);
    ba_problem.parameters_.push_back(pt3D[1]);
    ba_problem.parameters_.push_back(pt3D[2]);
  }

  // Create residuals for each observation in the bundle adjustment problem. The
  // parameters for cameras and points are added automatically.
  ceres::Problem problem;
  for (int i = 0; i < ba_problem.num_observations(); ++i) {

    // Each Residual block takes a point and a camera as input and outputs a 2
    // dimensional residual. Internally, the cost function stores the observed
    // image location and compares the reprojection against the observation.
    ceres::CostFunction* cost_function =
        new ceres::AutoDiffCostFunction<pinhole_reprojectionError::ErrorFunc_Refine_Intrinsic_Motion_3DPoints, 2, 3, 6, 3>(
            new pinhole_reprojectionError::ErrorFunc_Refine_Intrinsic_Motion_3DPoints(
                & ba_problem.observations()[2 * i]));

    problem.AddResidualBlock(cost_function,
                             NULL, // squared loss
                             ba_problem.mutable_camera_intrinsic_for_observation(i),
                             ba_problem.mutable_camera_extrinsic_for_observation(i),
                             ba_problem.mutable_point_for_observation(i));
  }

  // Make Ceres automatically detect the bundle structure. Note that the
  // standard solver, SPARSE_NORMAL_CHOLESKY, also works fine but it is slower
  // for standard bundle adjustment problems.
  ceres::Solver::Options options;
  options.linear_solver_type = ceres::SPARSE_SCHUR;
  if (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::SUITE_SPARSE))
    options.sparse_linear_algebra_library_type = ceres::SUITE_SPARSE;
  else
    if (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::CX_SPARSE))
      options.sparse_linear_algebra_library_type = ceres::CX_SPARSE;
    else
    {
      // No sparse backend for Ceres.
      // Use dense solving
      options.linear_solver_type = ceres::DENSE_SCHUR;
    }
  options.minimizer_progress_to_stdout = false;
  options.logging_type = ceres::SILENT;

  ceres::Solver::Summary summary;
  ceres::Solve(options, &problem, &summary);

  double dResidual_before = std::sqrt( summary.initial_cost / (ba_problem.num_observations_*2.));
  double dResidual_after = std::sqrt( summary.final_cost / (ba_problem.num_observations_*2.));

  std::cout << std::endl
    << "Bundle Adjustment of cameras [R|t], shared [f] and Structure : \n"
    << " Initial RMSE : " << dResidual_before << "\n"
    << " Final RMSE : " << dResidual_after << "\n"
    << "Initial focal : " << focal << "\n"
    << "Initial ppx : " << ppx << "\n"
    << "Initial ppy : " << ppy << std::endl;

  // If no error, get back refined parameters
  if (summary.IsSolutionUsable())
  {
    // Get back 3D points
    size_t cpt = 0;
    for (std::vector<Vec3>::iterator iter = vec_3DPoints.begin();
      iter != vec_3DPoints.end(); ++iter, ++cpt)
    {
      const double * pt = ba_problem.mutable_points() + cpt*3;
      Vec3 & pt3D = *iter;
      pt3D = Vec3(pt[0], pt[1], pt[2]);
    }
    // Get back camera
    for (cpt = 0; cpt < nCameraMotion; ++cpt)
    {
      const double * cam = ba_problem.mutable_cameras_extrinsic() + cpt*6;
      Mat3 R;
      // angle axis to rotation matrix
      ceres::AngleAxisToRotationMatrix(cam, R.data());
      Vec3 t(cam[3], cam[4], cam[5]);

      // Update the camera
      PinholeCamera & sCam = vec_cam[cpt];
      Mat3 K = sCam._K;
      double * intrinsics = ba_problem.mutable_cameras_intrinsic();
      std::cout << "Refined focal[" << cpt << "]: " << intrinsics[pinhole_reprojectionError::OFFSET_FOCAL_LENGTH] << std::endl;
      std::cout << "Refined ppx[" << cpt << "]: " << intrinsics[pinhole_reprojectionError::OFFSET_PRINCIPAL_POINT_X] << std::endl;
      std::cout << "Refined ppy[" << cpt << "]: " << intrinsics[pinhole_reprojectionError::OFFSET_PRINCIPAL_POINT_Y] << std::endl;
      K(0,0) = K(1,1) = focal;
      sCam = PinholeCamera(K, R, t);
    }
  }
}

/// Perform a Bundle Adjustment: Refine the cameras [R|t]
///  and common intrinsics [focal,ppx,ppy,k1,k2,k3] and the structure
void do_bundle_adjustment_common_intrinsics_brown_pinhole(
  PinholeCamera & camL,
  PinholeCamera & camR,
  const Mat & xL,
  const Mat & xR,
  const std::vector<size_t> & vec_inliers,
  std::vector<Vec3> & vec_3DPoints)
{
  int nCameraMotion = 2;
  int nCameraIntrinsic = 1;
  int n3Dpoints = vec_inliers.size();

  // Setup a BA problem
  BA_Problem_data_camMotionAndIntrinsic<6, 6> ba_problem;

  // Configure the size of the problem
  ba_problem.num_cameras_ = nCameraMotion;
  ba_problem.num_intrinsics_ = nCameraIntrinsic;
  ba_problem.num_points_ = n3Dpoints;
  ba_problem.num_observations_ = nCameraMotion * n3Dpoints;

  ba_problem.point_index_.reserve(ba_problem.num_observations_);
  ba_problem.camera_index_extrinsic.reserve(ba_problem.num_observations_);
  ba_problem.camera_index_intrinsic.reserve(ba_problem.num_observations_);
  ba_problem.observations_.reserve(2 * ba_problem.num_observations_);

  ba_problem.num_parameters_ =
    6 * ba_problem.num_cameras_ + 6 * ba_problem.num_intrinsics_ + 3 * ba_problem.num_points_;
  ba_problem.parameters_.reserve(ba_problem.num_parameters_);

  // Fill it with data (For each 3D point setup the tracks : the 2D visibility)
  // The two camera share the same intrinsic
  PinholeCamera vec_cam[2] = {camL, camR};
  for (int i = 0; i < n3Dpoints; ++i) {
    // Collect the image of point i in each frame (xL, xR).
    const Vec2 & xL_ = xL.col(vec_inliers[i]);
    const Vec2 & xR_ = xR.col(vec_inliers[i]);

    // Left 2D observations
    ba_problem.camera_index_extrinsic.push_back(0);
    ba_problem.camera_index_intrinsic.push_back(0);
    ba_problem.point_index_.push_back(i);
    ba_problem.observations_.push_back(xL_(0));
    ba_problem.observations_.push_back(xL_(1));

    // Right 2D observations
    ba_problem.camera_index_extrinsic.push_back(1);
    ba_problem.camera_index_intrinsic.push_back(0); // same intrinsic group
    ba_problem.point_index_.push_back(i);
    ba_problem.observations_.push_back(xR_(0));
    ba_problem.observations_.push_back(xR_(1));
  }

  // Add camera extrinsics [R,t]
  for (int j = 0; j < nCameraMotion; ++j) {
    // Rotation matrix to angle axis
    std::vector<double> angleAxis(3);
    ceres::RotationMatrixToAngleAxis((const double*)vec_cam[j]._R.data(), &angleAxis[0]);
    // translation
    Vec3 t = vec_cam[j]._t;
    ba_problem.parameters_.push_back(angleAxis[0]);
    ba_problem.parameters_.push_back(angleAxis[1]);
    ba_problem.parameters_.push_back(angleAxis[2]);
    ba_problem.parameters_.push_back(t[0]);
    ba_problem.parameters_.push_back(t[1]);
    ba_problem.parameters_.push_back(t[2]);
  }

  // Add camera intrinsics (focal, ppx, ppy, k1, k2, k3)
  {
    double focal = (vec_cam[0]._K(0,0) + vec_cam[0]._K(1,1)
      + vec_cam[1]._K(1,1) + vec_cam[1]._K(0,0)) / 4.0;
    double ppx = (vec_cam[0]._K(0,2) + vec_cam[1]._K(0,2)) / 2.0;
    double ppy = (vec_cam[0]._K(1,2) + vec_cam[1]._K(1,2)) / 2.0;
    double k1 = 0.0, k2 = 0.0, k3 = 0.0;

    // Setup intrinsics in the ba_problem
    ba_problem.parameters_.push_back(focal);
    ba_problem.parameters_.push_back(ppx);
    ba_problem.parameters_.push_back(ppy);
    ba_problem.parameters_.push_back(k1);
    ba_problem.parameters_.push_back(k2);
    ba_problem.parameters_.push_back(k3);
  }


  // Add 3D points coordinates parameters
  for (int i = 0; i < n3Dpoints; ++i) {
    Vec3 pt3D = vec_3DPoints[i];
    double * ptr3D = ba_problem.mutable_points() + i * 3;
    ba_problem.parameters_.push_back(pt3D[0]);
    ba_problem.parameters_.push_back(pt3D[1]);
    ba_problem.parameters_.push_back(pt3D[2]);
  }

  // Create residuals for each observation in the bundle adjustment problem. The
  // parameters for cameras and points are added automatically.
  ceres::Problem problem;
  for (int i = 0; i < ba_problem.num_observations(); ++i) {

    // Each Residual block takes a point and a camera as input and outputs a 2
    // dimensional residual. Internally, the cost function stores the observed
    // image location and compares the reprojection against the observation.
    ceres::CostFunction* cost_function =
        new ceres::AutoDiffCostFunction<pinhole_brown_reprojectionError::ErrorFunc_Refine_Camera_3DPoints, 2, 6, 6, 3>(
            new pinhole_brown_reprojectionError::ErrorFunc_Refine_Camera_3DPoints(
                & ba_problem.observations()[2 * i + 0]));

    problem.AddResidualBlock(cost_function,
                             NULL, // squared loss
                             ba_problem.mutable_camera_intrinsic_for_observation(i),
                             ba_problem.mutable_camera_extrinsic_for_observation(i),
                             ba_problem.mutable_point_for_observation(i));
  }

  // Make Ceres automatically detect the bundle structure. Note that the
  // standard solver, SPARSE_NORMAL_CHOLESKY, also works fine but it is slower
  // for standard bundle adjustment problems.
  ceres::Solver::Options options;
  options.linear_solver_type = ceres::SPARSE_SCHUR;
  if (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::SUITE_SPARSE))
    options.sparse_linear_algebra_library_type = ceres::SUITE_SPARSE;
  else
    if (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::CX_SPARSE))
      options.sparse_linear_algebra_library_type = ceres::CX_SPARSE;
    else
    {
      // No sparse backend for Ceres.
      // Use dense solving
      options.linear_solver_type = ceres::DENSE_SCHUR;
    }
  options.minimizer_progress_to_stdout = false;
  options.logging_type = ceres::SILENT;

  ceres::Solver::Summary summary;
  ceres::Solve(options, &problem, &summary);
  //std::cout << summary.FullReport() << std::endl;

  double dResidual_before = std::sqrt( summary.initial_cost / (ba_problem.num_observations_*2.));
  double dResidual_after = std::sqrt( summary.final_cost / (ba_problem.num_observations_*2.));

  std::cout << std::endl
    << "Bundle Adjustment of struture [X], cameras extrinsics [R|t],"
    << " and shared intrinsics [f,ppx,ppy,k1,k2,k3]: \n"
    << " Initial RMSE : " << dResidual_before << "\n"
    << " Final RMSE : " << dResidual_after << "\n"
    << "Refined intrinsics : " << std::endl;

  // If no error, get back refined parameters
  if (summary.IsSolutionUsable())
  {
    // Get back 3D points
    size_t cpt = 0;
    for (std::vector<Vec3>::iterator iter = vec_3DPoints.begin();
      iter != vec_3DPoints.end(); ++iter, ++cpt)
    {
      const double * pt = ba_problem.mutable_points() + cpt*3;
      Vec3 & pt3D = *iter;
      pt3D = Vec3(pt[0], pt[1], pt[2]);
    }
    // Get back camera
    for (cpt = 0; cpt < nCameraMotion; ++cpt)
    {
      const double * cam = ba_problem.mutable_cameras_extrinsic() + cpt*6;
      Mat3 R;
      // angle axis to rotation matrix
      ceres::AngleAxisToRotationMatrix(cam, R.data());
      Vec3 t(cam[3], cam[4], cam[5]);

      using namespace pinhole_brown_reprojectionError;
      // Update the camera
      PinholeCamera & sCam = vec_cam[cpt];
      const double * camIntrinsics = ba_problem.mutable_cameras_intrinsic();
      std::cout << " for camera Idx=[" << cpt << "]: " << std::endl
        << "\t focal: " << camIntrinsics[OFFSET_FOCAL_LENGTH] << std::endl
        << "\t ppx: " << camIntrinsics[OFFSET_PRINCIPAL_POINT_X] << std::endl
        << "\t ppy: " << camIntrinsics[OFFSET_PRINCIPAL_POINT_Y] << std::endl
        << "\t k1: " << camIntrinsics[OFFSET_K1] << std::endl
        << "\t k2: " << camIntrinsics[OFFSET_K2] << std::endl
        << "\t k3: " << camIntrinsics[OFFSET_K3] << std::endl
        << "\t initial: focal: " << sCam._K(0,0) << ", ppx: " << sCam._K(0,2)
        << ", ppy: " << sCam._K(1,2) <<std::endl;
      Mat3 K = sCam._K;
      K(0,0) = K(1,1) = camIntrinsics[OFFSET_FOCAL_LENGTH];
      sCam._K(0,2) = camIntrinsics[OFFSET_PRINCIPAL_POINT_X];
      sCam._K(1,2) = camIntrinsics[OFFSET_PRINCIPAL_POINT_Y];
      sCam = PinholeCamera(K, R, t);
    }
  }
}
