
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include "openMVG/image/image.hpp"
#include "openMVG/features/features.hpp"
#include "openMVG/matching/regions_matcher.hpp"
#include "openMVG/matching/svg_matches.hpp"
#include "openMVG/multiview/essential.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansac.hpp"
#include "openMVG/multiview/conditioning.hpp"
#include "openMVG/multiview/triangulation.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansacKernelAdaptator.hpp"

#include "nonFree/sift/SIFT_describer.hpp"
#include "openMVG_Samples/multiview_robust_essential_spherical/spherical_cam.hpp"

#include "software/SfM/SfMPlyHelper.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "third_party/vectorGraphics/svgDrawer.hpp"

#include <string>
#include <iostream>

using namespace openMVG;
using namespace openMVG::image;
using namespace openMVG::matching;
using namespace openMVG::robust;
using namespace svg;
using namespace std;

int main(int argc, char **argv) {

  CmdLine cmd;

  string jpg_filenameL, jpg_filenameR;

  cmd.add( make_option('a', jpg_filenameL, "input_a") );
  cmd.add( make_option('b', jpg_filenameR, "input_b") );

  std::cout << "Compute the relative pose between two spherical image."
   << "\nUse an Acontrario robust estimation based on angular errors." << std::endl;

  try
  {
    if (argc == 1)
    {
      const std::string sInputDir = std::string(THIS_SOURCE_DIR);
      jpg_filenameL = sInputDir + "/SponzaLion000.jpg";
      jpg_filenameR = sInputDir + "/SponzaLion001.jpg";
    }
    else
    {
      cmd.process(argc, argv);
    }
  } catch(const std::string& s) {

  }

  Image<unsigned char> imageL;
  ReadImage(jpg_filenameL.c_str(), &imageL);

  Image<unsigned char> imageR;
  ReadImage(jpg_filenameR.c_str(), &imageR);

  //--
  // Detect regions thanks to an image_describer
  //--
  using namespace openMVG::features;
  std::unique_ptr<Image_describer> image_describer
    (new SIFT_Image_describer(SIFT_Image_describer::Params(-1)));
  std::map<IndexT, std::unique_ptr<features::Regions> > regions_perImage;
  image_describer->Describe(imageL, regions_perImage[0]);
  image_describer->Describe(imageR, regions_perImage[1]);

  const SIFT_Regions* regionsL = dynamic_cast<SIFT_Regions*>(regions_perImage.at(0).get());
  const SIFT_Regions* regionsR = dynamic_cast<SIFT_Regions*>(regions_perImage.at(1).get());

  const PointFeatures
    featsL = regions_perImage.at(0)->GetRegionsPositions(),
    featsR = regions_perImage.at(1)->GetRegionsPositions();

  std::cout << "Left image SIFT count: " << featsL.size() << std::endl;
  std::cout << "Right image SIFT count: "<< featsR.size() << std::endl;

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
      regionsL->GetRegionsPositions(),
      jpg_filenameR,
      {imageR.Width(), imageR.Height()},
      regionsR->GetRegionsPositions(),
      "02_features.svg"
    );
  }

  std::vector<IndMatch> vec_PutativeMatches;
  //-- Perform matching -> find Nearest neighbor, filtered with Distance ratio
  {
    // Find corresponding points
    matching::DistanceRatioMatch(
      0.8, matching::ANN_L2,
      *regions_perImage.at(0).get(),
      *regions_perImage.at(1).get(),
      vec_PutativeMatches);

    IndMatchDecorator<float> matchDeduplicator(vec_PutativeMatches, featsL, featsR);
    matchDeduplicator.getDeduplicated(vec_PutativeMatches);

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

  // Essential geometry filtering of putative matches
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

    //-- Convert planar to spherical coordinates
    Mat xL_spherical(3,vec_PutativeMatches.size()), xR_spherical(3,vec_PutativeMatches.size());
    spherical_cam::planarToSpherical(xL, imageL.Width(), imageL.Height(), xL_spherical);
    spherical_cam::planarToSpherical(xR, imageR.Width(), imageR.Height(), xR_spherical);

    //-- Essential matrix robust estimation from spherical bearing vectors
    {
      std::vector<size_t> vec_inliers;

      // Define the AContrario angular error adaptor
      using KernelType =
        openMVG::robust::ACKernelAdaptor_AngularRadianError<
          // Use the 8 point solver in order to estimate E
          openMVG::spherical_cam::EightPointRelativePoseSolver,
          openMVG::spherical_cam::AngularError,
          Mat3>;

      KernelType kernel(xL_spherical, xR_spherical);

      // Robust estimation of the Essential matrix and it's precision
      Mat3 E;
      const double precision = std::numeric_limits<double>::infinity();
      const std::pair<double,double> ACRansacOut =
        ACRANSAC(kernel, vec_inliers, 1024, &E, precision, true);
      const double & threshold = ACRansacOut.first;

      std::cout << "\n Angular threshold found: " << R2D(threshold) << "(Degree)"<<std::endl;
      std::cout << "\n #Putatives/#inliers : " << xL_spherical.cols() << "/" << vec_inliers.size() << "\n" << std::endl;

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
        "04_inliers.svg",
        bVertical
      );

      if (vec_inliers.size() > 120)
      {
        // If an essential matrix have been found
        // Extract R|t
        //  - From the 4 possible solutions extracted from E keep the best
        //  - (check cheirality of correspondence once triangulation is done)

        // Accumulator to find the best solution
        std::vector<size_t> f(4, 0);

        std::vector<Mat3> Es;  // Essential,
        std::vector<Mat3> Rs;  // Rotation matrix.
        std::vector<Vec3> ts;  // Translation matrix.

        Es.push_back(E);
        // Recover best rotation and translation from E.
        MotionFromEssential(E, &Rs, &ts);

        //-> Test the 4 solutions will all the point
        Mat34 P1;
        P_From_KRt(Mat3::Identity(), Mat3::Identity(), Vec3::Zero(), &P1);
        std::vector< std::vector<size_t> > vec_newInliers(4);
        std::vector< std::vector<Vec3> > vec_3D(4);

        for (int kk = 0; kk < 4; ++kk) {
          const Mat3 &R2 = Rs[kk];
          const Vec3 &t2 = ts[kk];
          Mat34 P2;
          P_From_KRt(Mat3::Identity(), R2, t2, &P2);

          //-- For each inlier:
          //   - triangulate
          //   - check chierality
          for (const auto & inlier_it : vec_inliers)
          {
            const Vec3 & x1_ = xL_spherical.col(inlier_it);
            const Vec3 & x2_ = xR_spherical.col(inlier_it);

            //Triangulate
            Vec3 X;
            TriangulateDLT(P1, x1_, P2, x2_, &X);

            //Check positivity of the depth (sign of the dot product)
            const Vec3 Mc = R2 * X + t2;
            if (x2_.dot(Mc) > 0 && x1_.dot(X) > 0)
            {
              ++f[kk];
              vec_newInliers[kk].push_back(inlier_it);
              vec_3D[kk].push_back(X);
            }
          }
        }
        std::cout << std::endl << "estimate_Rt_fromE" << std::endl;
        std::cout << " #points in front of both cameras for each solution: "
          << f[0] << " " << f[1] << " " << f[2] << " " << f[3] << std::endl;

        std::vector<size_t>::iterator iter = max_element(f.begin(), f.end());
        if(*iter != 0)  {
          const size_t index = std::distance(f.begin(),iter);
          if (f[index] < 120) {
            std::cout << "Found solution have too few 3D points." << std::endl;
          }
          else  {
            std::cout << "Export found 3D scene in current directory." << std::endl;
            vec_inliers.clear();
            vec_inliers = vec_newInliers[index];
            std::ostringstream os;
            os << "./" << "relativePose_Spherical"<< ".ply";
            plyHelper::exportToPly(vec_3D[index], os.str());
          }
        }
      }
    }
  }
  return EXIT_SUCCESS;
}
