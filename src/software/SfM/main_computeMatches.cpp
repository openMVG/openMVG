
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/image/image.hpp"
#include "openMVG/features/features.hpp"
#include "openMVG/matching/matcher_brute_force.hpp"
#include "openMVG/matching/matcher_kdtree_flann.hpp"
#include "openMVG/matching/indMatch_utils.hpp"

/// Generic Image Collection image matching
#include "software/SfM/ImageCollectionMatcher_AllInMemory.hpp"
#include "software/SfM/ImageCollectionGeometricFilter.hpp"
#include "software/SfM/ImageCollection_F_ACRobust.hpp"
#include "software/SfM/pairwiseAdjacencyDisplay.hpp"

/// Feature detector and descriptor interface
#include "patented/sift/SIFT.hpp"

#include "third_party/cmdLine/cmdLine.h"

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>

using namespace openMVG;
using namespace openMVG::matching;
using namespace openMVG::robust;
using namespace std;

int main(int argc, char **argv)
{
  CmdLine cmd;

  std::string sImaDirectory;
  std::string sImaExtension;
  std::string sOutDir = "";
  float fDistRatio = .6f;
  bool bOctMinus1 = false;

  cmd.add( make_option('i', sImaDirectory, "imadir") );
  cmd.add( make_option('e', sImaExtension, "ext") );
  cmd.add( make_option('o', sOutDir, "outdir") );
  cmd.add( make_option('r', fDistRatio, "distratio") );
  cmd.add( make_option('s', bOctMinus1, "octminus1") );

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch(const std::string& s) {
      std::cerr << "Usage: " << argv[0] << ' '
      << "[-i|--imadir path] "
      << "[-e|--ext extension '*.jpg' or '*.png'] "
      << "[-o|--outdir path] "
      << "[-r|--distratio 0.6] "
      << "[-s|--octminus1 0 or 1] "
      << std::endl;

      std::cerr << s << std::endl;
      return EXIT_FAILURE;
  }

  std::cout << " You called : " <<std::endl
            << "main_computeMatches" << std::endl
            << "--imadir " << sImaDirectory << std::endl
            << "--ext " << sImaExtension << std::endl
            << "--outdir " << sOutDir << std::endl
            << "--distratio " << fDistRatio << std::endl
            << "--octminus1 " << bOctMinus1 << std::endl;

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

  //---------------------------------------
  // a. List images
  //---------------------------------------

  std::vector<std::string> vec_fileNames;
  std::vector<std::pair<size_t, size_t> > vec_imagesSize;

  if (!stlplus::folder_exists(sImaDirectory)) {
    std::cerr << "It is an invalid input image directory" << std::endl;
    return EXIT_FAILURE;
  } else  {

    //---------------------------------------
    // Look for images in the given directory
    //---------------------------------------
    vec_fileNames = stlplus::folder_wildcard(sImaDirectory,
      sImaExtension, false, true);

    std::sort(vec_fileNames.begin(), vec_fileNames.end());

    std::ofstream file(stlplus::create_filespec(sOutDir, "/lists", ".txt").c_str());
    std::copy(vec_fileNames.begin(), vec_fileNames.end(),
              std::ostream_iterator<std::string>(file, "\n"));
    file.close();

    for (size_t i=0; i < vec_fileNames.size(); ++i)  {
      vec_fileNames[i] = stlplus::create_filespec(
        stlplus::folder_append_separator(sImaDirectory), vec_fileNames[i]);
    }
    // DEBUG INFO
    std::cout << std::endl << "IMAGE(S) :" << std::endl;
    copy(vec_fileNames.begin(), vec_fileNames.end(), ostream_iterator<string>(cout, "\n"));

    if (vec_fileNames.empty())
    {
      std::cout << "\n No images in the provided directory.";
      return EXIT_FAILURE;
    }
  }


  //---------------------------------------
  // b. Compute features and descriptor
  //    - extract sift features and descriptor
  //    - if keypoints already computed, re-load them
  //    - else save features and descriptors on disk
  //---------------------------------------

  typedef Descriptor<float, 128> DescriptorT;
  typedef SIOPointFeature FeatureT;
  typedef std::vector<FeatureT> FeatsT;
  typedef vector<DescriptorT > DescsT;
  typedef KeypointSet<FeatsT, DescsT > KeypointSetT;

  // extract SIFT features
  {
    std::cout << "\n\nEXTRACT SIFT FEATURES" << std::endl;
    vec_imagesSize.resize(vec_fileNames.size());
    Image<RGBColor> imageRGB;
    Image<unsigned char> imageGray;

    C_Progress_display my_progress_bar( vec_fileNames.size() );
    for(size_t i=0; i < vec_fileNames.size(); ++i)  {
      KeypointSetT kpSet;

      std::string sFeat = stlplus::create_filespec(sOutDir,
        stlplus::basename_part(vec_fileNames[i]), "feat");
      std::string sDesc = stlplus::create_filespec(sOutDir,
        stlplus::basename_part(vec_fileNames[i]), "desc");

      //Test if descriptor and feature was already computed
      if (stlplus::file_exists(sFeat) && stlplus::file_exists(sDesc)) {

        if (ReadImage(vec_fileNames[i].c_str(), &imageRGB)) {
          vec_imagesSize[i] = make_pair(imageRGB.Width(), imageRGB.Height());
        }
        else {
          ReadImage(vec_fileNames[i].c_str(), &imageGray);
          vec_imagesSize[i] = make_pair(imageGray.Width(), imageGray.Height());
        }
      }
      else  { //Not already computed, so compute and save

        if (ReadImage(vec_fileNames[i].c_str(), &imageRGB)) {
          Rgb2Gray(imageRGB, &imageGray);
        }
        else{
          ReadImage(vec_fileNames[i].c_str(), &imageGray);
        }
        // Compute features and descriptors and export them to file
        SIFTDetector(imageGray,  kpSet.features(), kpSet.descriptors(), bOctMinus1);
        kpSet.saveToBinFile(sFeat, sDesc);
        vec_imagesSize[i] = make_pair(imageGray.Width(), imageRGB.Height());
      }
      ++my_progress_bar;
    }
  }



  //---------------------------------------
  // c. Compute putatives descriptor matches
  //    - L2 descriptor matching
  //    - Keep correspondences only if NearestNeighbor ratio is ok
  //---------------------------------------
  IndexedMatchPerPair map_PutativesMatches;
  // Define the matcher and the used metric (Squared L2)
  // ANN matcher could be defined as follow:
  typedef flann::L2<DescriptorT::bin_type> MetricT;
  typedef ArrayMatcher_Kdtree_Flann<DescriptorT::bin_type, MetricT> MatcherT;
  // Brute force matcher is defined as following:
  //typedef L2_Vectorized<DescriptorT::bin_type> MetricT;
  //typedef ArrayMatcherBruteForce<DescriptorT::bin_type, MetricT> MatcherT;

  // If the matches already exists, reload them
  if (stlplus::file_exists(sOutDir + "/matches.putative.txt"))
  {
    PairedIndMatchImport(sOutDir + "/matches.putative.txt", map_PutativesMatches);
    std::cout << std::endl << "PUTATIVE MATCHES -- PREVIOUS RESULTS LOADED" << std::endl;
  }
  else // Compute the putatives matches
  {
    ImageCollectionMatcher_AllInMemory<KeypointSetT, MatcherT> collectionMatcher(fDistRatio);
    if (collectionMatcher.loadData(vec_fileNames, sOutDir))
    {
      std::cout  << std::endl << "PUTATIVE MATCHES" << std::endl;
      collectionMatcher.Match(vec_fileNames, map_PutativesMatches);
      //---------------------------------------
      //-- Export putative matches
      //---------------------------------------
      std::ofstream file (std::string(sOutDir + "/matches.putative.txt").c_str());
      if (file.is_open())
        PairedIndMatchToStream(map_PutativesMatches, file);
      file.close();
    }
  }
  //-- export putative matches Adjacency matrix
  PairWiseMatchingToAdjacencyMatrixSVG(vec_fileNames.size(),
    map_PutativesMatches,
    stlplus::create_filespec(sOutDir, "PutativeAdjacencyMatrix", "svg"));


  //---------------------------------------
  // d. Geometric filtering of putatives matches
  //    - AContrario Estimation of the Fundamental matrix
  //    - Use a upper bound for the plausible F matrix
  //      acontrario estimated threshold
  //---------------------------------------
  IndexedMatchPerPair map_GeometricMatches_F;

  GeometricFilter_FMatrix_AC geomFilter_F_AC(4.0);
  ImageCollectionGeometricFilter<FeatureT, GeometricFilter_FMatrix_AC> collectionGeomFilter;
  if (collectionGeomFilter.loadData(vec_fileNames, sOutDir))
  {
    std::cout << std::endl << " - GEOMETRIC FILTERING - " << std::endl;
    collectionGeomFilter.Filter(
      geomFilter_F_AC,
      map_PutativesMatches,
      map_GeometricMatches_F,
      vec_imagesSize);

    //---------------------------------------
    //-- Export geometric filtered matches
    //---------------------------------------
    std::ofstream file (string(sOutDir + "/matches.f.txt").c_str());
    if (file.is_open())
      PairedIndMatchToStream(map_GeometricMatches_F, file);
    file.close();

    //-- export Adjacency matrix
    std::cout << "\n Export Adjacency Matrix of the pairwise's Epipolar matches"
      << std::endl;
    PairWiseMatchingToAdjacencyMatrixSVG(vec_fileNames.size(),
      map_GeometricMatches_F,
      stlplus::create_filespec(sOutDir, "EpipolarAdjacencyMatrix", "svg"));
  }
  return EXIT_SUCCESS;
}


