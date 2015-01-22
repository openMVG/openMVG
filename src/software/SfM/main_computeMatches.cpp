
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include "openMVG/image/image.hpp"
#include "openMVG/features/features.hpp"

/// Generic Image Collection image matching
#include "openMVG/matching_image_collection/Matcher_AllInMemory.hpp"
#include "openMVG/matching_image_collection/GeometricFilter.hpp"
#include "openMVG/matching_image_collection/F_ACRobust.hpp"
#include "openMVG/matching_image_collection/E_ACRobust.hpp"
#include "openMVG/matching_image_collection/H_ACRobust.hpp"
#include "software/SfM/pairwiseAdjacencyDisplay.hpp"
#include "software/SfM/SfMIOHelper.hpp"
#include "openMVG/matching/matcher_brute_force.hpp"
#include "openMVG/matching/matcher_kdtree_flann.hpp"
#include "openMVG/matching/indMatch_utils.hpp"

/// Feature detector and descriptor interface
#include "nonFree/sift/SIFT.hpp"

#include "openMVG/system/timer.hpp"

#include "third_party/cmdLine/cmdLine.h"

#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "third_party/progress/progress.hpp"

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>

using namespace openMVG;
using namespace openMVG::matching;
using namespace openMVG::robust;
using namespace std;

enum eGeometricModel
{
  FUNDAMENTAL_MATRIX = 0,
  ESSENTIAL_MATRIX   = 1,
  HOMOGRAPHY_MATRIX  = 2
};

enum ePairMode
{
  PAIR_EXHAUSTIVE = 0,
  PAIR_CONTIGUOUS = 1,
  PAIR_FROM_FILE  = 2,
  PAIR_PER_IMAGE = 3
};

bool readFile(const std::string& file, 
        std::vector<std::string>& vec_imageNames)
{
    std::ifstream in(file.c_str());
    if(!in.is_open())  {
      std::cerr << "Impossible to read file : " << file << std::endl;
      return false;
    }
    std::string imageName;
    while(std::getline(in, imageName))
        vec_imageNames.push_back(imageName);

    in.close();
    return true;
}

void prepareData(const std::vector<std::string>& vec_allImageNames,
        std::vector<std::string>& vec_imageNamesA,
        std::vector<std::string>& vec_imageNamesB,
        std::map<std::string, size_t>& map_imageNameToIndex,
        std::map<std::string, std::string>& map_imageNameToBaseName,
        std::map<std::string, std::vector<size_t> > & map_imageNameToComplementaryIndexes)
{
    // Compute indexes for list B
    for(std::vector<std::string>::const_iterator imageNameIt = vec_imageNamesB.begin(); imageNameIt != vec_imageNamesB.end(); ++imageNameIt)
    {
        // Compute index
        std::vector<std::string>::const_iterator imageName = find(vec_allImageNames.begin(), vec_allImageNames.end(), *imageNameIt);
        if(imageName == vec_allImageNames.end())
        {
            std::cerr << "\nImage \""<< *imageName << "\n not found." << std::endl;
            continue;
        }
        map_imageNameToIndex[*imageName] = std::distance<std::vector<std::string>::const_iterator>(vec_allImageNames.begin(), imageName);
    } 
    
    // Set data for list A
    for(std::vector<std::string>::const_iterator imageNameIt = vec_imageNamesA.begin(); imageNameIt != vec_imageNamesA.end(); ++imageNameIt)
    {
        // Compute index
        std::vector<std::string>::const_iterator imageName = find(vec_allImageNames.begin(), vec_allImageNames.end(), *imageNameIt);
        if(imageName == vec_allImageNames.end())
        {
            std::cerr << "\nImage \""<< *imageName << "\n not found." << std::endl;
            continue;
        }
        map_imageNameToIndex[*imageName] = std::distance<std::vector<std::string>::const_iterator>(vec_allImageNames.begin(), imageName);
        
        // Get basename
        map_imageNameToBaseName[*imageName] = stlplus::basename_part(*imageName) + ".";       
        
        // Compute complementary indexes
        std::vector<size_t> complementaryIndexes;
        for(std::vector<std::string>::const_iterator indexIt = vec_imageNamesB.begin(); indexIt != vec_imageNamesB.end(); ++indexIt)
        {
            size_t index = map_imageNameToIndex[*indexIt];
            if(map_imageNameToIndex[*imageNameIt] < index)
                complementaryIndexes.push_back(index);
        }
        map_imageNameToComplementaryIndexes[*imageNameIt] = complementaryIndexes;
        complementaryIndexes.clear();
         
    }
     
    assert(map_imageNameToBaseName.size() == map_imageNameToComplementaryIndexes.size());
}

int main(int argc, char **argv)
{
  CmdLine cmd;

  std::string sImageNameFileA;
  std::string sImageNameFileB;
  std::string sOutDir = "";
  std::string sGeometricModel = "f";
  float fDistRatio = .6f;
  int iMatchingVideoMode = -1;
  std::string sPredefinedPairList = "";

  cmd.add( make_option('o', sOutDir, "outdir") );
  cmd.add( make_option('a', sImageNameFileA, "imageNameFileA") );
  cmd.add( make_option('b', sImageNameFileB, "imageNameFileB") );
  cmd.add( make_option('r', fDistRatio, "distratio") );
  cmd.add( make_option('g', sGeometricModel, "geometricModel") );
  cmd.add( make_option('v', iMatchingVideoMode, "videoModeMatching") );
  cmd.add( make_option('l', sPredefinedPairList, "pairList") );

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch(const std::string& s) {
      std::cerr << "Usage: " << argv[0] << '\n'
      << "[-o|--outdir path] \n"
      << "\n[Optional]\n"
      << "[-f1|--imageNameFile1 filePath] \n"
      << "[-f2|--imageNameFile2 filePath] \n"
      << "[-r|--distratio 0.6] \n"
      << "[-g|--geometricModel f, e or h] \n"
      << "[-v|--videoModeMatching 2 -> X] \n"
      << "\t sequence matching with an overlap of X images\n"
      << "[-l]--pairList file"
      << std::endl;

      std::cerr << s << std::endl;
      return EXIT_FAILURE;
  }

  std::cout << " You called : " <<std::endl
            << argv[0] << std::endl
            << "--imageNameFileA " << sImageNameFileA << std::endl
            << "--imageNameFileB" << sImageNameFileB << std::endl
            << "--outdir " << sOutDir << std::endl
            << "--distratio " << fDistRatio << std::endl
            << "--geometricModel " << sGeometricModel << std::endl
            << "--videoModeMatching " << iMatchingVideoMode << std::endl;

  // TODO : check if ok with iMatchingVideoMode
  ePairMode ePairmode = (iMatchingVideoMode == -1 ) ? PAIR_PER_IMAGE : PAIR_CONTIGUOUS;

  if (sPredefinedPairList.length()) {
    std::cout << "--pairList " << sPredefinedPairList << std::endl;
    ePairmode = PAIR_FROM_FILE;
    if (iMatchingVideoMode>0) {
      std::cerr << "\nIncompatible options: --videoModeMatching and --pairList" << std::endl;
      return EXIT_FAILURE;
    }
  }

  if (sOutDir.empty())  {
    std::cerr << "\nIt is an invalid output directory" << std::endl;
    return EXIT_FAILURE;
  }

  // -----------------------------
  // a. List images
  // b. Compute putative descriptor matches
  // c. Geometric filtering of putative matches
  // -----------------------------

  // Create output dir
  if (!stlplus::folder_exists(sOutDir))
    stlplus::folder_create(sOutDir);

  //---------------------------------------
  // a. List images
  //---------------------------------------
  std::string sListsFile = stlplus::create_filespec(sOutDir, "lists.txt" );
  if (!stlplus::is_file(sListsFile)) {
    std::cerr << std::endl
      << "The input file \""<< sListsFile << "\" is missing" << std::endl;
    return false;
  }

  std::vector<openMVG::SfMIO::CameraInfo> vec_camImageName;
  std::vector<openMVG::SfMIO::IntrinsicCameraInfo> vec_focalGroup;
  if (!openMVG::SfMIO::loadImageList( vec_camImageName,
                                      vec_focalGroup,
                                      sListsFile) )
  {
    std::cerr << "\nEmpty or invalid image list." << std::endl;
    return false;
  }

  //-- Two alias to ease access to image filenames and image sizes
  std::vector<std::string> vec_allImageNames;
  std::vector<std::pair<size_t, size_t> > vec_imagesSize;
  for ( std::vector<openMVG::SfMIO::CameraInfo>::const_iterator
    iter_camInfo = vec_camImageName.begin();
    iter_camInfo != vec_camImageName.end();
    iter_camInfo++ )
  {
    vec_imagesSize.push_back( std::make_pair( vec_focalGroup[iter_camInfo->m_intrinsicId].m_w,
                                              vec_focalGroup[iter_camInfo->m_intrinsicId].m_h ) );
    vec_allImageNames.push_back(iter_camInfo->m_sImageName);
  }
    
  // Get image names for list 1
  std::vector<std::string> vec_imageNamesA;
  if(!sImageNameFileA.empty())
    readFile(sImageNameFileA, vec_imageNamesA);
  else
    vec_imageNamesA = vec_allImageNames;
  // Get image names for list 2
  std::vector<std::string> vec_imageNamesB;
  if(!sImageNameFileB.empty())
    readFile(sImageNameFileB, vec_imageNamesB);
  else
    vec_imageNamesB = vec_allImageNames;
  
  // TODO : remove indexes (use imageName as key)
  std::map<std::string, std::string> imageNameToBaseName;
  std::map<std::string, size_t> imageNameToIndex;
  std::map<std::string, std::vector<size_t> > imageNameToComplementaryIndexes;
  prepareData(vec_allImageNames, vec_imageNamesA, vec_imageNamesB, imageNameToIndex, imageNameToBaseName, imageNameToComplementaryIndexes);
  
  //---------------------------------------
  // b. Compute putative descriptor matches
  //    - L2 descriptor matching
  //    - Keep correspondences only if NearestNeighbor ratio is ok
  //---------------------------------------
  typedef Descriptor<unsigned char, 128> DescriptorT;
  typedef SIOPointFeature FeatureT;
  typedef std::vector<FeatureT> FeatsT;
  typedef vector<DescriptorT > DescsT;
  typedef KeypointSet<FeatsT, DescsT > KeypointSetT;
  PairWiseMatches map_PutativesMatches;
  // Define the matcher and the used metric (Squared L2)
  // ANN matcher could be defined as follow:
  typedef flann::L2<DescriptorT::bin_type> MetricT;
  typedef ArrayMatcher_Kdtree_Flann<DescriptorT::bin_type, MetricT> MatcherT;
  // Brute force matcher can be defined as following:
  //typedef L2_Vectorized<DescriptorT::bin_type> MetricT;
  //typedef ArrayMatcherBruteForce<DescriptorT::bin_type, MetricT> MatcherT;

  std::cout << std::endl << " - PUTATIVE MATCHES - " << std::endl;

  // TODO : add new matches to file (do no erase it))
  std::map<std::string, PairWiseMatches> imageNameToPairWiseMatches;
  for(std::vector<std::string>::const_iterator it = vec_imageNamesA.begin(); it != vec_imageNamesA.end(); ++it)
  {
    map_PutativesMatches.clear();
    const std::string& imageName = *it;
    std::string putativeFile = "/" + imageNameToBaseName[imageName] + "matches.putative.txt";
    if (stlplus::file_exists(sOutDir + putativeFile))
    {
      // Import files data in map_PutativesMatches
      PairedIndMatchImport(sOutDir + putativeFile, map_PutativesMatches);
      imageNameToPairWiseMatches[imageName] = map_PutativesMatches;
      std::cout << "\t PREVIOUS RESULTS LOADED" << std::endl;
    }
    else // Compute the putative matches
    {
      std::cout << "Use: ";
      switch (ePairmode)
      {
        case PAIR_EXHAUSTIVE: std::cout << "exhaustive pairwise matching" << std::endl; break;
        case PAIR_CONTIGUOUS: std::cout << "sequence pairwise matching" << std::endl; break;
        case PAIR_FROM_FILE:  std::cout << "user defined pairwise matching" << std::endl; break;
        case PAIR_PER_IMAGE:  std::cout << "pairwise matching per image" << std::endl; break;
      }    

      Timer timer;
      Matcher_AllInMemory<KeypointSetT, MatcherT> collectionMatcher(fDistRatio);
      // TODO : don't load all data
      if (collectionMatcher.loadData(vec_allImageNames, sOutDir)) // Load .feat and .desc for all images
      {
        // Get pair to match according the matching mode:
        PairsT pairs;
        switch (ePairmode)
        {
          case PAIR_EXHAUSTIVE: pairs = exhaustivePairs(vec_allImageNames.size()); break;
          case PAIR_CONTIGUOUS: pairs = contiguousWithOverlap(vec_allImageNames.size(), iMatchingVideoMode); break;
          case PAIR_PER_IMAGE: pairs = imagePairs(imageNameToIndex[imageName], imageNameToComplementaryIndexes[imageName]); break;
          case PAIR_FROM_FILE:
            if(!loadPairs(vec_allImageNames.size(), sPredefinedPairList, pairs))
            {
                return EXIT_FAILURE;
            };
            break;
        }
        // Photometric matching of putative pairs
        collectionMatcher.Match(vec_allImageNames, pairs, map_PutativesMatches);
        imageNameToPairWiseMatches[imageName] = map_PutativesMatches;
        //---------------------------------------
        //-- Export putative matches
        //---------------------------------------
        std::ofstream file (std::string(sOutDir + putativeFile).c_str());
        if (file.is_open())
          PairedIndMatchToStream(map_PutativesMatches, file);
        file.close();
      }
      std::cout << "Task done in (s): " << timer.elapsed() << std::endl;
    }
  }
  assert(imageNameToBaseName.size() == imageNameToPairWiseMatches.size());
  
  // TODO - On file per image ? 
  //-- export putative matches Adjacency matrix
//  PairWiseMatchingToAdjacencyMatrixSVG(vec_fileNames.size(),
//    map_PutativesMatches,
//    stlplus::create_filespec(sOutDir, "PutativeAdjacencyMatrix", "svg"));


  //---------------------------------------
  // c. Geometric filtering of putative matches
  //    - AContrario Estimation of the desired geometric model
  //    - Use an upper bound for the a contrario estimated threshold
  //---------------------------------------
  PairWiseMatches map_GeometricMatches;
   ImageCollectionGeometricFilter<FeatureT> collectionGeomFilter;
   const double maxResidualError = 4.0;
    // TODO : don't load data
  if (collectionGeomFilter.loadData(vec_allImageNames, sOutDir))
  {
    for(std::map<std::string, PairWiseMatches>::const_iterator it = imageNameToPairWiseMatches.begin(); it != imageNameToPairWiseMatches.end(); ++it)
    {
        Timer timer;
        const std::string& imageName = it->first;
        map_GeometricMatches.clear();
        std::cout << std::endl << " - GEOMETRIC FILTERING - " << std::endl;

        eGeometricModel eGeometricModelToCompute = FUNDAMENTAL_MATRIX;
        std::string sGeometricMatchesFilename = "";
        switch(sGeometricModel[0])
        {
          case 'f': case 'F':
            eGeometricModelToCompute = FUNDAMENTAL_MATRIX;
            sGeometricMatchesFilename = imageNameToBaseName[imageName] + "matches.f.txt";
          break;
          case 'e': case 'E':
            eGeometricModelToCompute = ESSENTIAL_MATRIX;
            sGeometricMatchesFilename = imageNameToBaseName[imageName] + "matches.e.txt";
          break;
          case 'h': case 'H':
            eGeometricModelToCompute = HOMOGRAPHY_MATRIX;
            sGeometricMatchesFilename = imageNameToBaseName[imageName] +"matches.h.txt";
          break;
          default:
            std::cerr << "Unknown geometric model" << std::endl;
            return EXIT_FAILURE;
        }

        switch (eGeometricModelToCompute)
        {
          case FUNDAMENTAL_MATRIX:
          {
            collectionGeomFilter.Filter(
               GeometricFilter_FMatrix_AC(maxResidualError),
               imageNameToPairWiseMatches[imageName],
               map_GeometricMatches,
               vec_imagesSize);
          }
          break;
          case ESSENTIAL_MATRIX:
          {
            // Build the intrinsic parameter map for each image
            std::map<size_t, Mat3> map_K;
            size_t cpt = 0;
            for ( std::vector<openMVG::SfMIO::CameraInfo>::const_iterator
              iter_camInfo = vec_camImageName.begin();
              iter_camInfo != vec_camImageName.end();
              ++iter_camInfo, ++cpt )
            {
              if (vec_focalGroup[iter_camInfo->m_intrinsicId].m_bKnownIntrinsic)
                map_K[cpt] = vec_focalGroup[iter_camInfo->m_intrinsicId].m_K;
            }

            collectionGeomFilter.Filter(
              GeometricFilter_EMatrix_AC(map_K, maxResidualError),
              imageNameToPairWiseMatches[imageName],
              map_GeometricMatches,
              vec_imagesSize);

            //-- Perform an additional check to remove pairs with poor overlap
            std::vector<PairWiseMatches::key_type> vec_toRemove;
            for (PairWiseMatches::const_iterator iterMap = map_GeometricMatches.begin();
              iterMap != map_GeometricMatches.end(); ++iterMap)
            {
              const size_t putativePhotometricCount = imageNameToPairWiseMatches[imageName].find(iterMap->first)->second.size();
              const size_t putativeGeometricCount = iterMap->second.size();
              const float ratio = putativeGeometricCount / (float)putativePhotometricCount;
              if (putativeGeometricCount < 50 || ratio < .3f)  {
                // the pair will be removed
                vec_toRemove.push_back(iterMap->first);
              }
            }
            //-- remove discarded pairs
            for (std::vector<PairWiseMatches::key_type>::const_iterator
              iter =  vec_toRemove.begin(); iter != vec_toRemove.end(); ++iter)
            {
              map_GeometricMatches.erase(*iter);
            }
          }
          break;
          case HOMOGRAPHY_MATRIX:
          {
            collectionGeomFilter.Filter(
              GeometricFilter_HMatrix_AC(maxResidualError),
              imageNameToPairWiseMatches[imageName],
              map_GeometricMatches,
              vec_imagesSize);
          }
          break;
        }
        
        //---------------------------------------
        //-- Export geometric filtered matches
        //---------------------------------------
        std::ofstream file (string(sOutDir + "/" + sGeometricMatchesFilename).c_str());
        if (file.is_open())
          PairedIndMatchToStream(map_GeometricMatches, file);
        file.close();
        std::cout << "Task done in (s): " << timer.elapsed() << std::endl;
    }

   
    //TODO : 
    //-- export Adjacency matrix
//    std::cout << "\n Export Adjacency Matrix of the pairwise's geometric matches"
//      << std::endl;
//    PairWiseMatchingToAdjacencyMatrixSVG(vec_fileNames.size(),
//      map_GeometricMatches,
//      stlplus::create_filespec(sOutDir, "GeometricAdjacencyMatrix", "svg"));
  }
  return EXIT_SUCCESS;
}
