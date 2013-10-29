
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include "openMVG/cameras/PinholeCamera.hpp"
#include "openMVG/image/image.hpp"
#include "openMVG/features/features.hpp"
#include "openMVG/matching/matcher_brute_force.hpp"
#include "openMVG/matching/indMatchDecoratorXY.hpp"
#include "openMVG/multiview/solver_essential_kernel.hpp"
#include "openMVG/multiview/projection.hpp"
#include "openMVG/multiview/triangulation.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansac.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansacKernelAdaptator.hpp"

#include "patented/sift/SIFT.hpp"
#include "openMVG_Samples/siftPutativeMatches/two_view_matches.hpp"

#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "third_party/vectorGraphics/svgDrawer.hpp"

#include <string>
#include <iostream>

using namespace openMVG;
using namespace openMVG::matching;
using namespace openMVG::robust;
using namespace svg;
using namespace std;

/// From the essential matrix test the 4 possible solutions
///  and return the best one. (point in front of the selected solution)
bool estimate_Rt_fromE(const Mat3 & K1, const Mat3 & K2,
  const Mat & x1, const Mat & x2,
  const Mat3 & E, const std::vector<size_t> & vec_inliers,
  Mat3 * R, Vec3 * t);

/// Read intrinsic K matrix from a file (ASCII)
/// F 0 ppx
/// 0 F ppy
/// 0 0 1
bool readIntrinsic(const std::string & fileName, Mat3 & K);

/// Export 3D point vector and camera position to PLY format
bool exportToPly(const std::vector<Vec3> & vec_points,
  const std::vector<Vec3> & vec_camPos,
  const std::string & sFileName);

/// Estimate the essential matrix from point matches and K matrices.
bool robustEssential(
  const Mat3 & K1, const Mat3 & K2,
  const Mat & x1, const Mat & x2,
  Mat3 * pE,
  std::vector<size_t> * pvec_inliers,
  const std::pair<size_t, size_t> & size_ima1,
  const std::pair<size_t, size_t> & size_ima2,
  double * errorMax,
  double * NFA,
  double precision = std::numeric_limits<double>::infinity());

/// Triangulate and export valid point as PLY (point in front of the cameras)
void triangulateAndSaveResult(
  const PinholeCamera & camL,
  const PinholeCamera & camR,
  const std::vector<size_t> & vec_inliers,
  const Mat & xL,
  const Mat & xR);

int main() {

  std::string sInputDir = stlplus::folder_up(string(THIS_SOURCE_DIR))
    + "/imageData/SceauxCastle/";
  Image<RGBColor> image;
  string jpg_filenameL = sInputDir + "100_7101.jpg";
  ReadImage(jpg_filenameL.c_str(), &image);

  Image<unsigned char> imageL;
  Rgb2Gray(image, &imageL);

  string jpg_filenameR = sInputDir + "100_7102.jpg";
  ReadImage(jpg_filenameR.c_str(), &image);

  Image<unsigned char> imageR;
  Rgb2Gray(image, &imageR);

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
    readIntrinsic(stlplus::create_filespec(sInputDir,"K","txt"), K);

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
        size_t idx = vec_inliers[i];
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

      //C. Triangulate and export inliers as a PLY scene
      triangulateAndSaveResult(
        camL, camR,
        vec_inliers,
        xL, xR);
    }
    else  {
      std::cout << "ACRANSAC was unable to estimate a rigid essential matrix"
        << std::endl;
    }
  }
  return EXIT_SUCCESS;
}

/// Estimate the best possible Rotation/Translation from E
bool estimate_Rt_fromE(const Mat3 & K1, const Mat3 & K2,
  const Mat & x1, const Mat & x2,
  const Mat3 & E, const std::vector<size_t> & vec_inliers,
  Mat3 * R, Vec3 * t)
{
  bool bOk = false;

  // Accumulator to find the best solution
  std::vector<size_t> f(4, 0);

  std::vector<Mat3> Es; // Essential,
  std::vector<Mat3> Rs;  // Rotation matrix.
  std::vector<Vec3> ts;  // Translation matrix.

  Es.push_back(E);
  // Recover best rotation and translation from E.
  MotionFromEssential(E, &Rs, &ts);

  //-> Test the 4 solutions will all the point
  assert(Rs.size() == 4);
  assert(ts.size() == 4);

  Mat34 P1, P2;
  Mat3 R1 = Mat3::Identity();
  Vec3 t1 = Vec3::Zero();
  P_From_KRt(K1, R1, t1, &P1);

  for (int i = 0; i < 4; ++i) {
    const Mat3 &R2 = Rs[i];
    const Vec3 &t2 = ts[i];
    P_From_KRt(K2, R2, t2, &P2);
    Vec3 X;

    for (size_t k = 0; k < vec_inliers.size(); ++k)
    {
      const Vec2 & x1_ = x1.col(vec_inliers[k]);
      const Vec2 & x2_ = x2.col(vec_inliers[k]);
      TriangulateDLT(P1, x1_, P2, x2_, &X);
      // Test if point is front to the two cameras.
      if (Depth(R1, t1, X) > 0 && Depth(R2, t2, X) > 0) {
          ++f[i];
      }
    }
  }
  // Check the solution :
  std::cout << "\t Number of points in front of both cameras:"
    << f[0] << " " << f[1] << " " << f[2] << " " << f[3] << std::endl;
  std::vector<size_t>::iterator iter = max_element(f.begin(), f.end());
  if(*iter != 0)
  {
    size_t index = std::distance(f.begin(), iter);
    (*R) = Rs[index];
    (*t) = ts[index];
    bOk = true;
  }
  else  {
    std::cout << std::endl << "/!\\There is no right solution,"
      <<" probably intermediate results are not correct or no points"
      <<" in front of both cameras" << std::endl;
    bOk = false;
  }
  return bOk;
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

/// Estimate the essential matrix from point matches and K matrices.
bool robustEssential(
  const Mat3 & K1, const Mat3 & K2,
  const Mat & x1, const Mat & x2,
  Mat3 * pE,
  std::vector<size_t> * pvec_inliers,
  const std::pair<size_t, size_t> & size_ima1,
  const std::pair<size_t, size_t> & size_ima2,
  double * errorMax,
  double * NFA,
  double precision)
{
  assert(pvec_inliers != NULL);
  assert(pE != NULL);

  // Use the 5 point solver to estimate E
  typedef openMVG::essential::kernel::FivePointKernel SolverType;
  // Define the AContrario adaptor
  typedef ACKernelAdaptorEssential<
      SolverType,
      openMVG::fundamental::kernel::EpipolarDistanceError,
      UnnormalizerT,
      Mat3>
      KernelType;

  KernelType kernel(x1, size_ima1.first, size_ima1.second,
                    x2, size_ima2.first, size_ima2.second, K1, K2);

  // Robustly estimation of the Essential matrix and it's precision
  std::pair<double,double> ACRansacOut = ACRANSAC(kernel, *pvec_inliers,
    4096, pE, precision, true);
  *errorMax = ACRansacOut.first;
  *NFA = ACRansacOut.second;

  return pvec_inliers->size() > 2.5 * SolverType::MINIMUM_SAMPLES;
}

/// Triangulate and export valid point as PLY (point in front of the cameras)
void triangulateAndSaveResult(
  const PinholeCamera & camL,
  const PinholeCamera & camR,
  const std::vector<size_t> & vec_inliers,
  const Mat & xL,
  const Mat & xR)
{
  std::vector<Vec3> vec_3DPoints;
  std::vector<double> vec_residuals;
  size_t nbPointWithNegativeDepth = 0;
  for (size_t k = 0; k < vec_inliers.size(); ++k) {
    const Vec2 & xL_ = xL.col(vec_inliers[k]);
    const Vec2 & xR_ = xR.col(vec_inliers[k]);

    Vec3 X = Vec3::Zero();
    TriangulateDLT(camL._P, xL_, camR._P, xR_, &X);

    // Compute residual:
    double dResidual = (camL.Residual(X, xL_) + camR.Residual(X, xR_))/2.0;
    vec_residuals.push_back(dResidual);
    if (camL.Depth(X) < 0 && camR.Depth(X) < 0) {
      ++nbPointWithNegativeDepth;
    }
    else  {
      vec_3DPoints.push_back(X);
    }
  }
  if (nbPointWithNegativeDepth > 0)
  {
    std::cout << nbPointWithNegativeDepth
      << " correspondence(s) with negative depth have been discarded."
      << std::endl;
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

  // Export as PLY (camera pos + 3Dpoints)
  std::vector<Vec3> vec_camPos;
  vec_camPos.push_back( camL._C );
  vec_camPos.push_back( camR._C );
  exportToPly(vec_3DPoints, vec_camPos, "EssentialGeometry.ply");
}
