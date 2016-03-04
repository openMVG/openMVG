
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/image/image.hpp"
#include "openMVG/features/akaze/AKAZE.hpp"

#include "openMVG/matching/matcher_brute_force.hpp"
#include "openMVG/matching/matching_filters.hpp"
#include "openMVG_Samples/siftPutativeMatches/two_view_matches.hpp"

#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "third_party/vectorGraphics/svgDrawer.hpp"

#include <string>
#include <iostream>

using namespace openMVG;
using namespace openMVG::matching;
using namespace svg;
using namespace std;

int main() {

  Image<RGBColor> image;
  const std::string jpg_filenameL =
  stlplus::folder_up(string(THIS_SOURCE_DIR))  + "/imageData/StanfordMobileVisualSearch/Ace_0.png";
  const std::string jpg_filenameR =
  stlplus::folder_up(string(THIS_SOURCE_DIR))  + "/imageData/StanfordMobileVisualSearch/Ace_1.png";

  Image<unsigned char> imageL, imageR;
  ReadImage(jpg_filenameL.c_str(), &imageL);
  ReadImage(jpg_filenameR.c_str(), &imageR);

  // Prepare vector to store detected feature and associated descriptor
  vector<SIOPointFeature> featsL, featsR;
  // Call AKAZE detector
  AKAZEConfig options;

#define AKAZE_MSURF 1
//#define AKAZE_MLBD 1
//#define AKAZE_LIOP 1

#if defined AKAZE_MSURF
  typedef Descriptor<float,64> Descriptor_T;
  vector<Descriptor_T > descsL, descsR;
  AKAZEDetector<Descriptor_T::bin_type, 64>(imageL, featsL, descsL, options);
  AKAZEDetector<Descriptor_T::bin_type, 64>(imageR, featsR, descsR, options);
#endif
#if defined AKAZE_MLBD
  typedef Descriptor<std::bitset<486>, 1> Descriptor_T;
  vector<Descriptor_T > descsL, descsR;
  AKAZEDetector<Descriptor_T::bin_type, 1>(imageL, featsL, descsL, options);
  AKAZEDetector<Descriptor_T::bin_type, 1>(imageR, featsR, descsR, options);
#endif
#if defined AKAZE_LIOP
  // LIOP
  typedef unsigned char descType;
  typedef Descriptor<descType,144> Descriptor_T;
  vector<Descriptor_T > descsL, descsR;
  AKAZEDetector<descType, 144>(imageL, featsL, descsL, options);
  AKAZEDetector<descType, 144>(imageR, featsR, descsR, options);
#endif


  // Show both images side by side
  {
    Image<unsigned char> concat;
    ConcatH(imageL, imageR, concat);
    string out_filename = "00_images.jpg";
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
    string out_filename = "01_features.jpg";
    WriteImage(out_filename.c_str(), concat);
  }

  //-- Perform matching -> find Nearest neighbor, filtered with Distance ratio
  std::vector<IndMatch> vec_PutativeMatches;
  {
    // Define the matcher (BruteForce) and the used metric (Squared L2)
    #if defined AKAZE_MSURF || AKAZE_LIOP
    typedef L2_Vectorized<Descriptor_T::bin_type> Metric;
    #endif
    #if defined AKAZE_MLBD // Binary based descriptor use the Hamming metric
    typedef HammingBitSet<Descriptor_T::bin_type> Metric;
    #endif
    typedef ArrayMatcherBruteForce<Descriptor_T::bin_type, Metric> MatcherT;

    // Distance ratio squared due to squared metric
    getPutativesMatches<Descriptor_T, MatcherT>(descsL, descsR, Square(0.8), vec_PutativeMatches);

    // Draw correspondences after Nearest Neighbor ratio filter
    svgDrawer svgStream( imageL.Width() + imageR.Width(), max(imageL.Height(), imageR.Height()));
    svgStream.drawImage(jpg_filenameL, imageL.Width(), imageL.Height());
    svgStream.drawImage(jpg_filenameR, imageR.Width(), imageR.Height(), imageL.Width());
    for (size_t i = 0; i < vec_PutativeMatches.size(); ++i) {
      //Get back linked feature, draw a circle and link them by a line
      const SIOPointFeature & L = featsL[vec_PutativeMatches[i].i_];
      const SIOPointFeature & R = featsR[vec_PutativeMatches[i].j_];
      svgStream.drawLine(L.x(), L.y(), R.x()+imageL.Width(), R.y(), svgStyle().stroke("green", 2.0));
      svgStream.drawCircle(L.x(), L.y(), L.scale(), svgStyle().stroke("yellow", 2.0));
      svgStream.drawCircle(R.x()+imageL.Width(), R.y(), R.scale(),svgStyle().stroke("yellow", 2.0));
    }
    string out_filename = "02_akaze_Matches.svg";
    ofstream svgFile( out_filename.c_str() );
    svgFile << svgStream.closeSvgFile().str();
    svgFile.close();
  }

  // Display some statistics
  std::cout
    << featsL.size() << " Features on image A" << std::endl
    << featsR.size() << " Features on image B" << std::endl
    << vec_PutativeMatches.size() << " matches after matching with Distance Ratio filter" << std::endl;

  return EXIT_SUCCESS;
}
