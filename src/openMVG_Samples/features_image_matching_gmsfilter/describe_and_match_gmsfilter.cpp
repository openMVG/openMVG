// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/image/image_io.hpp"
#include "openMVG/image/image_concat.hpp"
#include "openMVG/features/akaze/image_describer_akaze.hpp"
#include "openMVG/features/sift/SIFT_Anatomy_Image_Describer.hpp"
#include "openMVG/features/svg_features.hpp"
#include "openMVG/matching/regions_matcher.hpp"
#include "openMVG/matching/svg_matches.hpp"
#include "openMVG/robust_estimation/gms_filter.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <memory>
#include <string>

using namespace openMVG;
using namespace openMVG::image;
using namespace openMVG::robust;
using namespace std;

int main(int argc, char **argv) {

  // Add options to choose the desired Image_describer
  std::string sImage_describer_type = "SIFT";

  CmdLine cmd;
  cmd.add( make_option('t', sImage_describer_type, "type") );
  cmd.add( make_switch('d', "distance_ratio"));

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch(const std::string& s) {
      std::cerr << "Usage: " << argv[0] << '\n'
      << "\n[Optional]\n"
      << "[-t|--type\n"
      << "  (choose an image_describer interface):\n"
      << "   SIFT: SIFT keypoint & descriptor,\n"
      << "   AKAZE: AKAZE keypoint & floating point descriptor]\n"
      << "   AKAZE_MLDB: AKAZE keypoint & binary descriptor]\n"
      << "[-d|distance_ratio] Use distance ratio filter before GMS, else 1-1 matching is used."
      << std::endl;

      std::cerr << s << std::endl;
      return EXIT_FAILURE;
  }

  const string jpg_filenameL = stlplus::folder_up(string(THIS_SOURCE_DIR))
    + "/imageData/StanfordMobileVisualSearch/Ace_0.png";
  const string jpg_filenameR = stlplus::folder_up(string(THIS_SOURCE_DIR))
    + "/imageData/StanfordMobileVisualSearch/Ace_1.png";

  Image<unsigned char> imageL, imageR;
  ReadImage(jpg_filenameL.c_str(), &imageL);
  ReadImage(jpg_filenameR.c_str(), &imageR);
  assert(imageL.data() && imageR.data());

  // Call Keypoint extractor
  using namespace openMVG::features;
  std::unique_ptr<Image_describer> image_describer;
  if (sImage_describer_type == "SIFT")
    image_describer.reset(new SIFT_Anatomy_Image_describer(SIFT_Anatomy_Image_describer::Params()));
  else if (sImage_describer_type == "AKAZE")
    image_describer = AKAZE_Image_describer::create
      (AKAZE_Image_describer::Params(AKAZE::Params(), AKAZE_MSURF));
  else if (sImage_describer_type == "AKAZE_MLDB")
    image_describer = AKAZE_Image_describer::create
      (AKAZE_Image_describer::Params(AKAZE::Params(), AKAZE_MLDB));

  if (!image_describer)
  {
    std::cerr << "Invalid Image_describer type" << std::endl;
    return EXIT_FAILURE;
  }

  // GMS requires to have a lot of features
  image_describer->Set_configuration_preset(features::ULTRA_PRESET);

  //--
  // Detect regions thanks to the image_describer
  //--
  std::map<IndexT, std::unique_ptr<features::Regions> > regions_perImage;
  image_describer->Describe(imageL, regions_perImage[0]);
  image_describer->Describe(imageR, regions_perImage[1]);

  //--
  // Display used images & Features
  //--

  {
    //- Show images side by side
    Image<unsigned char> concat;
    ConcatH(imageL, imageR, concat);
    const string out_filename = "00_images.jpg";
    WriteImage(out_filename.c_str(), concat);
  }

  {
    //- Draw features on the images (side by side)
    Features2SVG
    (
      jpg_filenameL,
      {imageL.Width(), imageL.Height()},
      regions_perImage.at(0)->GetRegionsPositions(),
      jpg_filenameR,
      {imageR.Width(), imageR.Height()},
      regions_perImage.at(1)->GetRegionsPositions(),
      "01_features.svg"
    );
  }

  //--
  // Compute corresponding points
  //--
  //-- Perform matching -> find Nearest neighbor, filtered with Distance ratio
  matching::IndMatches vec_PutativeMatches;
  const bool distance_ratio_matching = cmd.used('d');
  if (distance_ratio_matching)
  {
    const float kDistanceRatio = 0.8f;
    matching::DistanceRatioMatch(
      kDistanceRatio,
      (sImage_describer_type == "AKAZE_MLDB") ? matching::BRUTE_FORCE_HAMMING
        : matching::BRUTE_FORCE_L2,
      *regions_perImage.at(0).get(),
      *regions_perImage.at(1).get(),
      vec_PutativeMatches);
  }
  else
  {
    matching::Match(
      (sImage_describer_type == "AKAZE_MLDB") ? matching::BRUTE_FORCE_HAMMING
        : matching::BRUTE_FORCE_L2,
      *regions_perImage.at(0).get(),
      *regions_perImage.at(1).get(),
      vec_PutativeMatches);
  }
  // Draw the putative photometric correspondences
  {
    const bool bVertical = true;
    Matches2SVG
    (
      jpg_filenameL,
      {imageL.Width(), imageL.Height()},
      regions_perImage.at(0)->GetRegionsPositions(),
      jpg_filenameR,
      {imageR.Width(), imageR.Height()},
      regions_perImage.at(1)->GetRegionsPositions(),
      vec_PutativeMatches,
      "02_Matches.svg",
      bVertical
    );
  }

  // Display some statistics
  std::cout << vec_PutativeMatches.size() << " #matches Found" << std::endl;

  // Apply the GMS filter
  {
    std::vector<Eigen::Vector2f> vec_points_left;
    {
      const auto & regions_pos = regions_perImage.at(0)->GetRegionsPositions();
      vec_points_left.reserve(regions_pos.size());
      for (const auto & it : regions_pos)
        vec_points_left.push_back({it.x(), it.y()});
    }
    std::vector<Eigen::Vector2f> vec_points_right;
    {
      const auto & regions_pos = regions_perImage.at(1)->GetRegionsPositions();
      vec_points_right.reserve(regions_pos.size());
      for (const auto & it : regions_pos)
        vec_points_right.push_back({it.x(), it.y()});
    }

    const int kGmsThreshold = 6;
    robust::GMSFilter gms(
      vec_points_left,  {imageL.Width(), imageL.Height()},
      vec_points_right, {imageR.Width(), imageR.Height()},
      vec_PutativeMatches,
      kGmsThreshold);
    const bool with_scale_invariance = true;
    const bool with_rotation_invariance = true;
    std::vector<bool> inlier_flags;
    const int nb_inliers = gms.GetInlierMask(inlier_flags,
                                             with_scale_invariance,
                                             with_rotation_invariance);

    std::cout
      << vec_points_left.size() << " #Features on image A" << std::endl
      << vec_points_right.size() << " #Features on image B" << std::endl
      << nb_inliers << " #matches kept by the GMS Filter" << std::endl;

    matching::IndMatches vec_gms_matches;
    for (int i = 0; i < static_cast<int>(inlier_flags.size()); ++i)
    {
      if (inlier_flags[i])
        vec_gms_matches.push_back(vec_PutativeMatches[i]);
    }
    // Draw the correspondences kept by the GMSFilter
    {
      const bool bVertical = true;
      Matches2SVG
      (
        jpg_filenameL,
        {imageL.Width(), imageL.Height()},
        regions_perImage.at(0)->GetRegionsPositions(),
        jpg_filenameR,
        {imageR.Width(), imageR.Height()},
        regions_perImage.at(1)->GetRegionsPositions(),
        vec_gms_matches,
        "03_GMSMatches.svg",
        bVertical
      );
    }
  }
  return EXIT_SUCCESS;
}
