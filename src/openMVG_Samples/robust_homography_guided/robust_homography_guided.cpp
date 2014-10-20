
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/image/image.hpp"
#include "openMVG/features/features.hpp"
#include "openMVG/matching/matcher_brute_force.hpp"
#include "openMVG/multiview/solver_homography_kernel.hpp"
#include "openMVG/multiview/conditioning.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansac.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansacKernelAdaptator.hpp"

#include "openMVG/robust_estimation/guided_matching.hpp"

#include "nonFree/sift/SIFT.hpp"
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

int main() {

  Image<RGBColor> image;
  string jpg_filenameL = stlplus::folder_up(string(THIS_SOURCE_DIR))
    + "/imageData/StanfordMobileVisualSearch/Ace_0.png";
  string jpg_filenameR = stlplus::folder_up(string(THIS_SOURCE_DIR))
    + "/imageData/StanfordMobileVisualSearch/Ace_1.png";

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
  bool bOctaveMinus1 = true;
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

  // Homography geometry filtering of putative matches
  {
    //A. get back interest point and send it to the robust estimation framework
    Mat xL(2, vec_PutativeMatches.size());
    Mat xR(2, vec_PutativeMatches.size());

    for (size_t k = 0; k < vec_PutativeMatches.size(); ++k)  {
      const SIOPointFeature & imaL = featsL[vec_PutativeMatches[k]._i];
      const SIOPointFeature & imaR = featsR[vec_PutativeMatches[k]._j];
      xL.col(k) = imaL.coords().cast<double>();
      xR.col(k) = imaR.coords().cast<double>();
    }

    //-- Homography robust estimation
    std::vector<size_t> vec_inliers;
    typedef ACKernelAdaptor<
      openMVG::homography::kernel::FourPointSolver,
      openMVG::homography::kernel::AsymmetricError,
      UnnormalizerI,
      Mat3>
      KernelType;

    KernelType kernel(
      xL, imageL.Width(), imageL.Height(),
      xR, imageR.Width(), imageR.Height(),
      false); // configure as point to point error model.

    Mat3 H;
    std::pair<double,double> ACRansacOut = ACRANSAC(kernel, vec_inliers, 1024, &H,
      std::numeric_limits<double>::infinity(),
      true);
    const double & thresholdH = ACRansacOut.first;

    // Check the homography support some point to be considered as valid
    if (vec_inliers.size() > KernelType::MINIMUM_SAMPLES *2.5) {

      std::cout << "\nFound a homography under the confidence threshold of: "
        << thresholdH << " pixels\n\twith: " << vec_inliers.size() << " inliers"
        << " from: " << vec_PutativeMatches.size()
         << " putatives correspondences"
        << std::endl;

      //Show homography validated point and compute residuals
      std::vector<double> vec_residuals(vec_inliers.size(), 0.0);
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
        // residual computation
        vec_residuals[i] = std::sqrt(KernelType::ErrorT::Error(H,
                                       LL.coords().cast<double>(),
                                       RR.coords().cast<double>()));
      }
      string out_filename = "04_ACRansacHomography.svg";
      ofstream svgFile( out_filename.c_str() );
      svgFile << svgStream.closeSvgFile().str();
      svgFile.close();

      // Display some statistics of reprojection errors
      float dMin, dMax, dMean, dMedian;
      minMaxMeanMedian<float>(vec_residuals.begin(), vec_residuals.end(),
                            dMin, dMax, dMean, dMedian);

      std::cout << std::endl
        << "Homography matrix estimation, residuals statistics:" << "\n"
        << "\t-- Residual min:\t" << dMin << std::endl
        << "\t-- Residual median:\t" << dMedian << std::endl
        << "\t-- Residual max:\t "  << dMax << std::endl
        << "\t-- Residual mean:\t " << dMean << std::endl;

      // Perform GUIDED MATCHING
      // Use the computed model with all possible point couples
      //  and keep the one that have the error under the AC-RANSAC precision
      //  value.

      Mat xL, xR;
      PointsToMat(featsL, xL);
      PointsToMat(featsR, xR);
      std::vector<IndMatch> vec_corresponding_index;
      GuidedMatching<Mat3, openMVG::homography::kernel::AsymmetricError>(
        H, xL, xR, Square(thresholdH), vec_corresponding_index);

      std::cout << "\nGuided homography matching found "
        << vec_corresponding_index.size() << " correspondences."
        << std::endl;

      // Merge AC Ransac and H-Guided matches
      std::set<IndMatch> set_matches(vec_corresponding_index.begin(),
        vec_corresponding_index.end());
      for ( size_t i = 0; i < vec_inliers.size(); ++i)  {
        set_matches.insert(vec_PutativeMatches[vec_inliers[i]]);
      }
      std::cout << "\nGuided homography + AC Ransac results matching found "
        << set_matches.size() << " correspondences."
        << std::endl;
      // Update the corresponding index:
      vec_corresponding_index.clear();
      vec_corresponding_index.assign(set_matches.begin(), set_matches.end());

      {
        //Show homography validated correspondences and compute residuals
        std::vector<double> vec_residuals(vec_corresponding_index.size(), 0.0);
        svgDrawer svgStream( imageL.Width() + imageR.Width(), max(imageL.Height(), imageR.Height()));
        svgStream.drawImage(jpg_filenameL, imageL.Width(), imageL.Height());
        svgStream.drawImage(jpg_filenameR, imageR.Width(), imageR.Height(), imageL.Width());
        for ( size_t i = 0; i < vec_corresponding_index.size(); ++i)  {

          const SIOPointFeature & LL = featsL[vec_corresponding_index[i]._i];
          const SIOPointFeature & RR = featsR[vec_corresponding_index[i]._j];
          const Vec2f L = LL.coords();
          const Vec2f R = RR.coords();
          svgStream.drawLine(L.x(), L.y(), R.x()+imageL.Width(), R.y(), svgStyle().stroke("green", 2.0));
          svgStream.drawCircle(L.x(), L.y(), LL.scale(), svgStyle().stroke("yellow", 2.0));
          svgStream.drawCircle(R.x()+imageL.Width(), R.y(), RR.scale(),svgStyle().stroke("yellow", 2.0));
          // residual computation
          vec_residuals[i] = std::sqrt(KernelType::ErrorT::Error(H,
                                         LL.coords().cast<double>(),
                                         RR.coords().cast<double>()));
        }
        string out_filename = "04_ACRansacHomography_guided.svg";
        ofstream svgFile( out_filename.c_str() );
        svgFile << svgStream.closeSvgFile().str();
        svgFile.close();

        // Display some statistics of reprojection errors
        float dMin, dMax, dMean, dMedian;
        minMaxMeanMedian<float>(vec_residuals.begin(), vec_residuals.end(),
                              dMin, dMax, dMean, dMedian);

        std::cout << std::endl
          << "Homography matrix estimation, residuals statistics:" << "\n"
          << "\t-- Residual min:\t" << dMin << std::endl
          << "\t-- Residual median:\t" << dMedian << std::endl
          << "\t-- Residual max:\t "  << dMax << std::endl
          << "\t-- Residual mean:\t " << dMean << std::endl;
      }
    }
    else  {
      std::cout << "ACRANSAC was unable to estimate a rigid homography"
        << std::endl;
    }
  }
  return EXIT_SUCCESS;
}
