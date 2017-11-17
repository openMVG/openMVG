// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/features/feature.hpp"
#include "openMVG/features/sift/SIFT_Anatomy_Image_Describer.hpp"
#include "openMVG/features/svg_features.hpp"
#include "openMVG/image/image_io.hpp"
#include "openMVG/image/image_concat.hpp"
#include "openMVG/numeric/eigen_alias_definition.hpp"
#include "openMVG/matching/regions_matcher.hpp"
#include "openMVG/matching/svg_matches.hpp"
#include "openMVG/multiview/solver_homography_kernel.hpp"
#include "openMVG/robust_estimation/guided_matching.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansac.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansacKernelAdaptator.hpp"
#include "openMVG/types.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "third_party/vectorGraphics/svgDrawer.hpp"

#include <cstdlib>
#include <iostream>
#include <string>
#include <utility>

using namespace openMVG;
using namespace openMVG::image;
using namespace openMVG::matching;
using namespace openMVG::robust;
using namespace svg;
using namespace std;

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
  using namespace openMVG::features;
  std::unique_ptr<Image_describer> image_describer
    (new SIFT_Anatomy_Image_describer(SIFT_Anatomy_Image_describer::Params(-1)));
  std::map<IndexT, std::unique_ptr<features::Regions>> regions_perImage;
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

  // Homography geometry filtering of putative matches
  {
    //A. get back interest point and send it to the robust estimation framework
    Mat xL(2, vec_PutativeMatches.size());
    Mat xR(2, vec_PutativeMatches.size());

    for (size_t k = 0; k < vec_PutativeMatches.size(); ++k)  {
      const PointFeature & imaL = featsL[vec_PutativeMatches[k].i_];
      const PointFeature & imaR = featsR[vec_PutativeMatches[k].j_];
      xL.col(k) = imaL.coords().cast<double>();
      xR.col(k) = imaR.coords().cast<double>();
    }

    //-- Homography robust estimation
    std::vector<uint32_t> vec_inliers;
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

    Mat3 H;
    const std::pair<double,double> ACRansacOut = ACRANSAC(kernel, vec_inliers, 1024, &H,
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
      const bool bVertical = true;
      InlierMatches2SVG
      (
        jpg_filenameL,
        {imageL.Width(), imageL.Height()},
        regionsL->GetRegionsPositions(),
        jpg_filenameR,
        {imageR.Width(), imageR.Height()},
        regionsR->GetRegionsPositions(),
        vec_PutativeMatches,
        vec_inliers,
        "04_ACRansacHomography.svg",
        bVertical
      );
      std::vector<double> vec_residuals(vec_inliers.size(), 0.0);
      for ( size_t i = 0; i < vec_inliers.size(); ++i)  {
        const SIOPointFeature & LL = regionsL->Features()[vec_PutativeMatches[vec_inliers[i]].i_];
        const SIOPointFeature & RR = regionsR->Features()[vec_PutativeMatches[vec_inliers[i]].j_];
        // residual computation
        vec_residuals[i] = std::sqrt(KernelType::ErrorT::Error(H,
                                       LL.coords().cast<double>(),
                                       RR.coords().cast<double>()));
      }

      // Display some statistics of reprojection errors
      float dMin, dMax, dMean, dMedian;
      minMaxMeanMedian<float>(vec_residuals.cbegin(), vec_residuals.cend(),
                            dMin, dMax, dMean, dMedian);

      std::cout << std::endl
        << "Homography matrix estimation, residuals statistics:" << "\n"
        << "\t-- Residual min:\t" << dMin << std::endl
        << "\t-- Residual median:\t" << dMedian << std::endl
        << "\t-- Residual max:\t "  << dMax << std::endl
        << "\t-- Residual mean:\t " << dMean << std::endl;

      // --
      // Perform GUIDED MATCHING
      // --
      // Use the computed model to check valid correspondences
      // a. by considering only the geometric error,
      // b. by considering geometric error and descriptor distance ratio.
      std::vector<IndMatches> vec_corresponding_indexes(2);

      Mat xL, xR;
      PointsToMat(featsL, xL);
      PointsToMat(featsR, xR);

      //a. by considering only the geometric error

      geometry_aware::GuidedMatching<Mat3, openMVG::homography::kernel::AsymmetricError>(
        H, xL, xR, Square(thresholdH), vec_corresponding_indexes[0]);
      std::cout << "\nGuided homography matching (geometric error) found "
        << vec_corresponding_indexes[0].size() << " correspondences."
        << std::endl;

      // b. by considering geometric error and descriptor distance ratio
      geometry_aware::GuidedMatching
        <Mat3, openMVG::homography::kernel::AsymmetricError>(
        H,
        nullptr, *regions_perImage.at(0), // nullptr since no Intrinsic is defined
        nullptr, *regions_perImage.at(1), // nullptr since no Intrinsic is defined
        Square(thresholdH), Square(0.8),
        vec_corresponding_indexes[1]);

      std::cout << "\nGuided homography matching "
        << "(geometric + descriptor distance ratio) found "
        << vec_corresponding_indexes[1].size() << " correspondences."
        << std::endl;

      for (size_t idx = 0; idx < 2; ++idx)
      {
        const std::vector<IndMatch> & vec_corresponding_index = vec_corresponding_indexes[idx];
        //Show homography validated correspondences
        svgDrawer svgStream( imageL.Width() + imageR.Width(), max(imageL.Height(), imageR.Height()));
        svgStream.drawImage(jpg_filenameL, imageL.Width(), imageL.Height());
        svgStream.drawImage(jpg_filenameR, imageR.Width(), imageR.Height(), imageL.Width());
        for ( size_t i = 0; i < vec_corresponding_index.size(); ++i)  {

          const SIOPointFeature & LL = regionsL->Features()[vec_corresponding_index[i].i_];
          const SIOPointFeature & RR = regionsR->Features()[vec_corresponding_index[i].j_];
          const Vec2f L = LL.coords();
          const Vec2f R = RR.coords();
          svgStream.drawLine(L.x(), L.y(), R.x()+imageL.Width(), R.y(), svgStyle().stroke("green", 2.0));
          svgStream.drawCircle(L.x(), L.y(), LL.scale(), svgStyle().stroke("yellow", 2.0));
          svgStream.drawCircle(R.x()+imageL.Width(), R.y(), RR.scale(),svgStyle().stroke("yellow", 2.0));
        }
        const string out_filename =
          (idx == 0) ? "04_ACRansacHomography_guided_geom.svg"
            : "04_ACRansacHomography_guided_geom_distratio.svg";
        ofstream svgFile( out_filename.c_str() );
        svgFile << svgStream.closeSvgFile().str();
        svgFile.close();
      }
    }
    else  {
      std::cout << "ACRANSAC was unable to estimate a rigid homography"
        << std::endl;
    }
  }
  return EXIT_SUCCESS;
}
