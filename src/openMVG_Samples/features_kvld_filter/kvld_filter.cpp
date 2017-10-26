// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 openMVG authors.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/features/sift/SIFT_Anatomy_Image_Describer.hpp"
#include "openMVG/features/svg_features.hpp"
#include "openMVG/image/image_io.hpp"
#include "openMVG/image/image_concat.hpp"
#include "openMVG/matching/kvld/kvld.h"
#include "openMVG/matching/kvld/kvld_draw.h"
#include "openMVG/matching/regions_matcher.hpp"
#include "openMVG/matching/svg_matches.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "third_party/vectorGraphics/svgDrawer.hpp"

#include <cstdlib>
#include <iostream>
#include <string>

using namespace openMVG;
using namespace openMVG::image;
using namespace openMVG::matching;
using namespace svg;
using namespace std;

int main(int argc, char **argv) {
  CmdLine cmd;

  std::string sImg1 = stlplus::folder_up(string(THIS_SOURCE_DIR))
    + "/imageData/StanfordMobileVisualSearch/Ace_0.png";
  std::string sImg2 = stlplus::folder_up(string(THIS_SOURCE_DIR))
    + "/imageData/StanfordMobileVisualSearch/Ace_1.png";
  std::string sOutDir = "./kvldOut";
  std::cout << sImg1 << std::endl << sImg2 << std::endl;
  cmd.add( make_option('i', sImg1, "img1") );
  cmd.add( make_option('j', sImg2, "img2") );
  cmd.add( make_option('o', sOutDir, "outdir") );

  if (argc > 1)
  {
    try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
    } catch (const std::string& s) {
        std::cerr << "Usage: " << argv[0] << ' '
        << "[-i|--img1 file] "
        << "[-j|--img2 file] "
        << "[-o|--outdir path] "
        << std::endl;

        std::cerr << s << std::endl;
        return EXIT_FAILURE;
    }
  }

  std::cout << " You called : " <<std::endl
            << argv[0] << std::endl
            << "--img1 " << sImg1 << std::endl
            << "--img2 " << sImg2 << std::endl
            << "--outdir " << sOutDir << std::endl;

  if (sOutDir.empty())  {
    std::cerr << "\nIt is an invalid output directory" << std::endl;
    return EXIT_FAILURE;
  }


  // -----------------------------
  // a. List images
  // b. Compute features and descriptor
  // c. Compute putatives descriptor matches
  // d. Geometric filtering of putatives matches
  // e. Export some statistics
  // -----------------------------

  // Create output dir
  if (!stlplus::folder_exists(sOutDir))
    stlplus::folder_create( sOutDir );

  const string jpg_filenameL = sImg1;
  const string jpg_filenameR = sImg2;

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
    string out_filename = "00_images.jpg";
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
      "01_features.svg"
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
      "02_Matches.svg",
      bVertical,
      std::max(std::max(imageL.Width(), imageL.Height()) / float(600), 2.0f)
    );
  }

  //K-VLD filter
  Image<float> imgA (imageL.GetMat().cast<float>());
  Image<float> imgB (imageR.GetMat().cast<float>());

  std::vector<Pair> matchesFiltered;
  std::vector<Pair> matchesPair;

  for (const auto & match_it : vec_PutativeMatches)
  {
    matchesPair.emplace_back(match_it.i_, match_it.j_);
  }
  std::vector<double> vec_score;

  //In order to illustrate the gvld(or vld)-consistant neighbors,
  // the following two parameters has been externalized as inputs of the function KVLD.
  openMVG::Mat E = openMVG::Mat::Ones(vec_PutativeMatches.size(), vec_PutativeMatches.size())*(-1);
  // gvld-consistancy matrix, intitialized to -1,  >0 consistancy value, -1=unknow, -2=false
  std::vector<bool> valid(vec_PutativeMatches.size(), true);// indices of match in the initial matches, if true at the end of KVLD, a match is kept.

  size_t it_num=0;
  KvldParameters kvldparameters; // initial parameters of KVLD
  while (it_num < 5 &&
          kvldparameters.inlierRate > KVLD(imgA, imgB, regionsL->Features(), regionsR->Features(),
          matchesPair, matchesFiltered, vec_score,E,valid,kvldparameters)) {
    kvldparameters.inlierRate /= 2;
    //std::cout<<"low inlier rate, re-select matches with new rate="<<kvldparameters.inlierRate<<std::endl;
    kvldparameters.K = 2;
    it_num++;
  }

  std::vector<IndMatch> vec_FilteredMatches;
  for (std::vector<Pair>::const_iterator i_matchFilter = matchesFiltered.begin();
      i_matchFilter != matchesFiltered.end(); ++i_matchFilter){
    vec_FilteredMatches.push_back(IndMatch(i_matchFilter->first, i_matchFilter->second));
  }

  //Print K-VLD consistent matches
  {
    svgDrawer svgStream(imageL.Width() + imageR.Width(), max(imageL.Height(), imageR.Height()));

    // ".svg"
    svgStream.drawImage(jpg_filenameL, imageL.Width(), imageL.Height());
    svgStream.drawImage(jpg_filenameR, imageR.Width(), imageR.Height(), imageL.Width());


    for (size_t it1=0; it1<matchesPair.size()-1;it1++){
      for (size_t it2=it1+1; it2<matchesPair.size();it2++){
         if (valid[it1] && valid[it2] && E(it1,it2)>=0){

          const PointFeature & l1 = featsL[matchesPair[it1].first];
          const PointFeature & r1 = featsR[matchesPair[it1].second];

          const PointFeature & l2 = featsL[matchesPair[it2].first];
          const PointFeature & r2 = featsR[matchesPair[it2].second];

          // Compute the width of the current VLD segment
          float L = (l1.coords() - l2.coords()).norm();
          float width = std::max(1.f, L / (dimension+1.f));

          // ".svg"
          svgStream.drawLine(l1.x(), l1.y(), l2.x(), l2.y(), svgStyle().stroke("yellow", width));
          svgStream.drawLine(r1.x() + imageL.Width(), r1.y(), r2.x() + imageL.Width(), r2.y(), svgStyle().stroke("yellow", width));

        }
      }
    }
    const string out_filename = stlplus::create_filespec(sOutDir, "03_KVLD_Matches.svg");
    ofstream svgFile( out_filename.c_str() );
    svgFile << svgStream.closeSvgFile().str();
    svgFile.close();
  }

  {
    //Print keypoints kept by K-VLD
    svgDrawer svgStream(imageL.Width() + imageR.Width(), max(imageL.Height(), imageR.Height()));

    // ".svg"
    svgStream.drawImage(jpg_filenameL, imageL.Width(), imageL.Height());
    svgStream.drawImage(jpg_filenameR, imageR.Width(), imageR.Height(), imageL.Width());

    for (size_t it=0; it<matchesPair.size();it++){
       if (valid[it]){

        const PointFeature & l = featsL[matchesPair[it].first];
        const PointFeature & r = featsR[matchesPair[it].second];

        // ".svg"
        svgStream.drawCircle(l.x(), l.y(), 10, svgStyle().stroke("yellow", 2.0));
        svgStream.drawCircle(r.x() + imageL.Width(), r.y(), 10, svgStyle().stroke("yellow", 2.0));
      }
    }
    const string out_filename = stlplus::create_filespec(sOutDir, "04_KVLD_Keypoints.svg");
    ofstream svgFile( out_filename.c_str() );
    svgFile << svgStream.closeSvgFile().str();
    svgFile.close();
  }

  Image <unsigned char> imageOutL = imageL;
  Image <unsigned char> imageOutR = imageR;

  getKVLDMask(
    &imageOutL, &imageOutR,
    regionsL->Features(), regionsR->Features(),
    matchesPair,
    valid,
    E);

  {
    const string out_filename = stlplus::create_filespec(sOutDir, "05_Left-K-VLD-MASK.jpg");
    WriteImage(out_filename.c_str(), imageOutL);
  }
  {
    const string out_filename = stlplus::create_filespec(sOutDir, "06_Right-K-VLD-MASK.jpg");
    WriteImage(out_filename.c_str(), imageOutR);
  }

  return EXIT_SUCCESS;
}
