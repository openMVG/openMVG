// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2016 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/image/image_io.hpp"
#include "openMVG/image/image_concat.hpp"
#include "openMVG/features/feature.hpp"
#include "openMVG/features/sift/SIFT_Anatomy_Image_Describer.hpp"
#include "openMVG/features/svg_features.hpp"
#include "openMVG/matching/regions_matcher.hpp"
#include "openMVG/matching/svg_matches.hpp"

// Robust estimation includes
//--
//-- Homography estimator
#include "openMVG/multiview/solver_homography_kernel.hpp"
//--
//- Fitting metric helper (measure the fitting error between the datum and a Model)
#include "openMVG/robust_estimation/score_evaluator.hpp"
//--
//- MaxConsensus
#include "openMVG/robust_estimation/robust_estimator_MaxConsensus.hpp"
//- Ransac
#include "openMVG/robust_estimation/robust_estimator_Ransac.hpp"
//- LMeds
#include "openMVG/robust_estimation/robust_estimator_LMeds.hpp"
//- ACRANSAC
#include "openMVG/robust_estimation/robust_estimator_ACRansac.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansacKernelAdaptator.hpp"

#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "third_party/vectorGraphics/svgDrawer.hpp"

#include <string>
#include <iostream>

using namespace openMVG;
using namespace openMVG::image;
using namespace openMVG::features;
using namespace openMVG::matching;
using namespace openMVG::robust;
using namespace svg;
using namespace std;

//- Helper function
// Display a visual report of the inliers
void display_info
(
  const std::string & jpg_filenameL,
  const Image<unsigned char> & imageL,
  const PointFeatures & featsL,
  const std::string & jpg_filenameR,
  const Image<unsigned char> & imageR,
  const PointFeatures & featsR,
  const IndMatches & vec_PutativeMatches,
  const Mat3 &H,
  const std::vector<uint32_t> & vec_inliers,
  const std::string & sMethod
);

int main() {

  Image<RGBColor> image;
  const string jpg_filenameL = stlplus::folder_up(string(THIS_SOURCE_DIR))
    + "/imageData/StanfordMobileVisualSearch/Ace_0.png";
  const string jpg_filenameR = stlplus::folder_up(string(THIS_SOURCE_DIR))
    + "/imageData/StanfordMobileVisualSearch/Ace_1.png";

  Image<unsigned char> imageL, imageR;
  ReadImage(jpg_filenameL.c_str(), &imageL);
  ReadImage(jpg_filenameR.c_str(), &imageR);

  //--
  // Detect regions thanks to an image_describer
  //--
  std::unique_ptr<Image_describer> image_describer
    (new SIFT_Anatomy_Image_describer(SIFT_Anatomy_Image_describer::Params(-1)));
  std::map<IndexT, std::unique_ptr<features::Regions> > regions_perImage;
  image_describer->Describe(imageL, regions_perImage[0]);
  image_describer->Describe(imageR, regions_perImage[1]);

  const SIFT_Regions* regionsL = dynamic_cast<SIFT_Regions*>(regions_perImage.at(0).get());
  const SIFT_Regions* regionsR = dynamic_cast<SIFT_Regions*>(regions_perImage.at(1).get());

  const PointFeatures
    featsL = regions_perImage.at(0)->GetRegionsPositions(),
    featsR = regions_perImage.at(1)->GetRegionsPositions();

  // Show both images side by side
  {
    Image<unsigned char> concat;
    ConcatH(imageL, imageR, concat);
    string out_filename = "01_concat.jpg";
    WriteImage(out_filename.c_str(), concat);
  }

  //- Draw features on the two image (side by side)
  {
    Features2SVG
    (
      jpg_filenameL,
      {imageL.Width(), imageL.Height()},
      regionsL->Features(),
      jpg_filenameR,
      {imageR.Width(), imageR.Height()},
      regionsR->Features(),
      "02_features.svg"
    );
  }

  std::vector<IndMatch> vec_PutativeMatches;
  //-- Perform matching -> find Nearest neighbor, filtered with Distance ratio
  {
    // Find corresponding points
    matching::DistanceRatioMatch(
      0.8, matching::BRUTE_FORCE_L2,
      *regions_perImage.at(0).get(),
      *regions_perImage.at(1).get(),
      vec_PutativeMatches);

    // Draw correspondences after Nearest Neighbor ratio filter
    const bool bVertical = true;
    Matches2SVG
    (
      jpg_filenameL,
      {imageL.Width(), imageL.Height()},
      regionsL->GetRegionsPositions(),
      jpg_filenameR,
      {imageR.Width(), imageR.Height()},
      regionsR->GetRegionsPositions(),
      vec_PutativeMatches,
      "03_Matches.svg",
      bVertical
    );
  }

  //---
  // Homography geometry filtering of putative matches
  // - Show how to use the robust_estimation framework with different robust_estimator methods.
  //---

  // First we list the SIFT photometric corresponding points to Mat arrays (The datum).
  Mat xL(2, vec_PutativeMatches.size());
  Mat xR(2, vec_PutativeMatches.size());

  for (size_t k = 0; k < vec_PutativeMatches.size(); ++k)
  {
    // For each correspondence, add the Right & Left feature point positions
    const PointFeature & imaL = featsL[vec_PutativeMatches[k].i_];
    const PointFeature & imaR = featsR[vec_PutativeMatches[k].j_];
    xL.col(k) = imaL.coords().cast<double>();
    xR.col(k) = imaR.coords().cast<double>();
  }

  // Then we use a robust_estimator to find if a model can be fitted in the defined datum

  //--
  //-- Max Consensus
  //- Return the Model that have the most of inliers
  //- Perform all the iterations (no early exit)
  {
    std::cout
      << "----------------------------------\n"
      << "MAXConsensus -- Robust estimation \n"
      << "----------------------------------\n";

    // Define the Homography Kernel
    using KernelType = homography::kernel::UnnormalizedKernel;
    KernelType kernel(xL, xR);

    // The Model type
    Mat3 H;
    // The inlier list
    std::vector<uint32_t> vec_inliers;
    H =
      MaxConsensus
      (
        kernel, // The Kernel (embed the correspondences, the Model Solver & the fitting metric.
        ScorerEvaluator<KernelType>(4.0), // Upper bound of the tolerated error for the inliers
        &vec_inliers, // Inlier list
        1024 // max_iteration count
      );
    display_info(
      jpg_filenameL, imageL, featsL,
      jpg_filenameR, imageR, featsR,
      vec_PutativeMatches,
      H, vec_inliers, "MAXConsensus");
  }

  //--
  //-- Ransac
  //- Return the Model that have the most of inliers
  //- Remaining iteration count can be decreased if meaningful models are found
  {
    std::cout
      << "----------------------------------\n"
      << "RANSAC -- Robust estimation \n"
      << "----------------------------------\n";

     // Define the Homography Kernel
    using KernelType = homography::kernel::UnnormalizedKernel;
    KernelType kernel(xL, xR);

    // The Model type
    Mat3 H;
    // The inlier list
    std::vector<uint32_t> vec_inliers;
    // Inlier count
    size_t inlier_count = 0;
    H =
      RANSAC
      (
        kernel, // The Kernel (embed the correspondences, the Model Solver & the fitting metric.
        ScorerEvaluator<KernelType>(4.0), // Upper bound of the tolerated error for the inliers
        &vec_inliers, // Inlier list
        &inlier_count // Inlier count
      );
    display_info(
      jpg_filenameL, imageL, featsL,
      jpg_filenameR, imageR, featsR,
      vec_PutativeMatches,
      H, vec_inliers, "RANSAC");
  }

  //--
  //-- LMeds
  //- Return the Model that have the most of inliers
  //- Minimize the median of the error (no input precision to provide)
  {
    std::cout
      << "----------------------------------\n"
      << "LMEDS -- Robust estimation \n"
      << "----------------------------------\n";

    // Define the Homography Kernel
    using KernelType = homography::kernel::UnnormalizedKernel;
    KernelType kernel(xL, xR);

    // The Model type
    Mat3 H;
    double lmeds_threshold = 0.0;
    LeastMedianOfSquares
      (
        kernel, // The Kernel (embed the correspondences, the Model Solver & the fitting metric.
        &H, // Update model
        &lmeds_threshold // Upper bound of the tolerated error for the inliers
      );
    // List the inliers
    std::vector<uint32_t> vec_inliers;
    for (int i = 0; i < xL.cols(); ++i)
    {
      const double residual_error = std::sqrt(kernel.Error(i, H));
      if (residual_error < lmeds_threshold)
        vec_inliers.push_back(i);
    }
    display_info(
      jpg_filenameL, imageL, featsL,
      jpg_filenameR, imageR, featsR,
      vec_PutativeMatches,
      H, vec_inliers, "LMEDS");
  }

  //--
  //-- AContrario Ransac
  //- Apriori inlier/outlier error threshold is only an option.
  //- The robust estimator finds the best meaningful model and its associated precision
  {
    std::cout
      << "----------------------------------\n"
      << "ACRansac -- Robust estimation \n"
      << "----------------------------------\n";

    // Define the AContrario Kernel adaptor for the Homography Model
    using KernelType =
      ACKernelAdaptor<
        openMVG::homography::kernel::FourPointSolver,
        openMVG::homography::kernel::AsymmetricError,
        UnnormalizerI,
        Mat3>;

    KernelType kernel(
      xL, imageL.Width(), imageL.Height(),
      xR, imageR.Width(), imageR.Height(),
      false); // configure as point to point error model.

    // The Model type
    Mat3 H;
    // The inlier list
    std::vector<uint32_t> vec_inliers;
    // Call the Robust Estimator on the KERNEL
    const std::pair<double,double> ACRansacOut = // Return the precision & the associated NFA
      ACRANSAC(
        kernel, // The Kernel (embed the correspondences, the Model Solver & the fitting metric.
        vec_inliers, // The inlier list
        1024, // Max iteration count (it can be less if a meaningful is found at the first iterations)
        &H, // Returned model
        std::numeric_limits<double>::infinity(), // No apriori Threshold
        false // Verbose to console
       );
    const double & thresholdH = ACRansacOut.first;
    std::cout << "AContrario Upperbound error estimation is: " << thresholdH << std::endl;
    display_info(
      jpg_filenameL, imageL, featsL,
      jpg_filenameR, imageR, featsR,
      vec_PutativeMatches,
      H, vec_inliers, "ACRansac");

  }
  return EXIT_SUCCESS;
}

void display_info
(
  const std::string & jpg_filenameL,
  const Image<unsigned char> & imageL,
  const PointFeatures & featsL,
  const std::string & jpg_filenameR,
  const Image<unsigned char> & imageR,
  const PointFeatures & featsR,
  const IndMatches & vec_PutativeMatches,
  const Mat3 &H,
  const std::vector<uint32_t> & vec_inliers,
  const std::string & sMethod
)
{
  std::cout
    << "\nFound a homography with: " << vec_inliers.size() << " inliers"
    << " from: " << vec_PutativeMatches.size()
    << " putatives correspondences"
    << std::endl;

  //Show homography validated point and compute residuals
  std::vector<double> vec_residuals(vec_inliers.size(), 0.0);
  svgDrawer svgStream( imageL.Width() + imageR.Width(), max(imageL.Height(), imageR.Height()));
  svgStream.drawImage(jpg_filenameL, imageL.Width(), imageL.Height());
  svgStream.drawImage(jpg_filenameR, imageR.Width(), imageR.Height(), imageL.Width());
  for ( size_t i = 0; i < vec_inliers.size(); ++i)  {
    const PointFeature & L = featsL[vec_PutativeMatches[vec_inliers[i]].i_];
    const PointFeature & R = featsR[vec_PutativeMatches[vec_inliers[i]].j_];
    svgStream.drawLine(L.x(), L.y(), R.x()+imageL.Width(), R.y(), svgStyle().stroke("green", 2.0));
    svgStream.drawCircle(L.x(), L.y(), 5, svgStyle().stroke("yellow", 2.0));
    svgStream.drawCircle(R.x()+imageL.Width(), R.y(), 5,svgStyle().stroke("yellow", 2.0));
    // residual computation
    using KernelType = homography::kernel::UnnormalizedKernel;
    vec_residuals[i] = std::sqrt(
        KernelType::ErrorT::Error(H,
          L.coords().cast<double>(),
          R.coords().cast<double>()));
  }
  std::ostringstream os;
  os << sMethod << "_robust_fitting.svg";
  ofstream svgFile( os.str().c_str() );
  svgFile << svgStream.closeSvgFile().str();
  svgFile.close();

  // Display some statistics of reprojection errors
  double dMin, dMax, dMean, dMedian;
  minMaxMeanMedian<double>(vec_residuals.begin(), vec_residuals.end(),
                        dMin, dMax, dMean, dMedian);

  std::cout << std::endl
    << "Homography matrix estimation, residuals statistics:" << "\n"
    << "\t-- Residual min:\t" << dMin << std::endl
    << "\t-- Residual median:\t" << dMedian << std::endl
    << "\t-- Residual max:\t "  << dMax << std::endl
    << "\t-- Residual mean:\t " << dMean << std::endl;
}
