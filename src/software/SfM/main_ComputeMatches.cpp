
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/pipelines/sfm_engine.hpp"
#include "openMVG/sfm/pipelines/sfm_features_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_regions_provider.hpp"

/// Generic Image Collection image matching
#include "openMVG/matching_image_collection/Matcher_Regions_AllInMemory.hpp"
#include "openMVG/matching_image_collection/Cascade_Hashing_Matcher_Regions_AllInMemory.hpp"
#include "openMVG/matching_image_collection/GeometricFilter.hpp"
#include "openMVG/matching_image_collection/F_ACRobust.hpp"
#include "openMVG/matching_image_collection/E_ACRobust.hpp"
#include "openMVG/matching_image_collection/H_ACRobust.hpp"
#include "openMVG/matching/pairwiseAdjacencyDisplay.hpp"
#include "openMVG/matching/indMatch_utils.hpp"
#include "openMVG/system/timer.hpp"

#include "openMVG/graph/graph.hpp"
#include "openMVG/stl/stl.hpp"
#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <cstdlib>
#include <fstream>

using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::matching;
using namespace openMVG::robust;
using namespace openMVG::sfm;
using namespace openMVG::matching_image_collection;
using namespace std;

enum EGeometricModel
{
  FUNDAMENTAL_MATRIX = 0,
  ESSENTIAL_MATRIX   = 1,
  HOMOGRAPHY_MATRIX  = 2
};

enum EPairMode
{
  PAIR_EXHAUSTIVE = 0,
  PAIR_CONTIGUOUS = 1,
  PAIR_FROM_FILE  = 2
};

/// Compute corresponding features between a series of views:
/// - Load view images description (regions: features & descriptors)
/// - Compute putative local feature matches (descriptors matching)
/// - Compute geometric coherent feature matches (robust model estimation from putative matches)
/// - Export computed data
int main(int argc, char **argv)
{
  CmdLine cmd;

  std::string sSfM_Data_Filename;
  std::string sMatchesDirectory = "";
  std::string sGeometricModel = "f";
  float fDistRatio = 0.8f;
  int iMatchingVideoMode = -1;
  std::string sPredefinedPairList = "";
  int rangeStart = -1;
  int rangeSize = 0;
  bool bUpRight = false;
  std::string sNearestMatchingMethod = "AUTO";
  bool bForce = false;
  bool bGuided_matching = false;
  int imax_iteration = 2048;
  bool matchFilePerImage = false;
  bool orderPairs = false;

  //required
  cmd.add( make_option('i', sSfM_Data_Filename, "input_file") );
  cmd.add( make_option('o', sMatchesDirectory, "out_dir") );
  // Options
  cmd.add( make_option('r', fDistRatio, "ratio") );
  cmd.add( make_option('g', sGeometricModel, "geometric_model") );
  cmd.add( make_option('v', iMatchingVideoMode, "video_mode_matching") );
  cmd.add( make_option('l', sPredefinedPairList, "pair_list") );
  cmd.add( make_option('s', rangeStart, "range_start") );
  cmd.add( make_option('r', rangeSize, "range_size") );
  cmd.add( make_option('n', sNearestMatchingMethod, "nearest_matching_method") );
  cmd.add( make_option('f', bForce, "force") );
  cmd.add( make_option('m', bGuided_matching, "guided_matching") );
  cmd.add( make_option('I', imax_iteration, "max_iteration") );
  cmd.add( make_option('x', matchFilePerImage, "match_file_per_image") );
  cmd.add( make_option('p', orderPairs, "order_pairs") );

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch(const std::string& s) {
      std::cerr << "Usage: " << argv[0] << '\n'
      << "[-i|--input_file] a SfM_Data file\n"
      << "[-o|--out_dir path] output path where computed matches are stored\n"
      << "\n[Optional]\n"
      << "[-f|--force] Force to recompute data]\n"
      << "[-r|--ratio] Distance ratio to discard non meaningful matches\n"
      << "   0.8: (default).\n"
      << "[-g|--geometric_model]\n"
      << "  (pairwise correspondences filtering thanks to robust model estimation):\n"
      << "   f: (default) fundamental matrix,\n"
      << "   e: essential matrix,\n"
      << "   h: homography matrix.\n"
      << "[-v|--video_mode_matching]\n"
      << "  (sequence matching with an overlap of X images)\n"
      << "   X: with match 0 with (1->X), ...]\n"
      << "   2: will match 0 with (1,2), 1 with (2,3), ...\n"
      << "   3: will match 0 with (1,2,3), 1 with (2,3,4), ...\n"
      << "[-l]--pair_list] filepath\n"
      << "  A file which contains the list of matches to perform.\n"
      << "[-s]--range_start] range image index start\n"
      << "  To compute only the matches for specified range.\n"
      << "  This allows to compute different matches on different computers in parallel.\n"
      << "[-r]--range_size] range size\n"
      << "  To compute only the matches for specified range.\n"
      << "  This allows to compute different matches on different computers in parallel.\n"
      << "[-n|--nearest_matching_method]\n"
      << "  AUTO: auto choice from regions type,\n"
      << "  For Scalar based regions descriptor:\n"
      << "    BRUTEFORCEL2: L2 BruteForce matching,\n"
      << "    ANNL2: L2 Approximate Nearest Neighbor matching,\n"
      << "    CASCADEHASHINGL2: L2 Cascade Hashing matching.\n"
      << "    FASTCASCADEHASHINGL2: (default)\n"
      << "      L2 Cascade Hashing with precomputed hashed regions\n"
      << "     (faster than CASCADEHASHINGL2 but use more memory).\n"
      << "  For Binary based descriptor:\n"
      << "    BRUTEFORCEHAMMING: BruteForce Hamming matching.\n"
      << "[-m|--guided_matching]\n"
      << "  use the found model to improve the pairwise correspondences.\n"
      << "[-x|--match_file_per_image]\n"
      << "  Save matches in a separate file per image\n"
      << "[-p|--order_pairs]\n"
      << "  Order the pairs"
      << std::endl;

      std::cerr << s << std::endl;
      return EXIT_FAILURE;
  }

  std::cout << " You called : " << "\n"
            << argv[0] << "\n"
            << "--input_file " << sSfM_Data_Filename << "\n"
            << "--out_dir " << sMatchesDirectory << "\n"
            << "Optional parameters:" << "\n"
            << "--force " << bForce << "\n"
            << "--ratio " << fDistRatio << "\n"
            << "--geometric_model " << sGeometricModel << "\n"
            << "--video_mode_matching " << iMatchingVideoMode << "\n"
            << "--pair_list " << sPredefinedPairList << "\n"
            << "--range_start " << rangeStart <<  "\n"
            << "--range_size " << rangeSize <<  "\n"
            << "--nearest_matching_method " << sNearestMatchingMethod << "\n"
            << "--guided_matching " << bGuided_matching << "\n"
            << "--order_pairs " << orderPairs << "\n"
            << "--match_file_per_image " << matchFilePerImage << std::endl;

  EPairMode ePairmode = (iMatchingVideoMode == -1 ) ? PAIR_EXHAUSTIVE : PAIR_CONTIGUOUS;

  if (sPredefinedPairList.length()) {
    ePairmode = PAIR_FROM_FILE;
    if (iMatchingVideoMode>0) {
      std::cerr << "\nIncompatible options: --videoModeMatching and --pairList" << std::endl;
      return EXIT_FAILURE;
    }
  }

  if (sMatchesDirectory.empty() || !stlplus::is_folder(sMatchesDirectory))  {
    std::cerr << "\nIt is an invalid output directory" << std::endl;
    return EXIT_FAILURE;
  }

  EGeometricModel eGeometricModelToCompute = FUNDAMENTAL_MATRIX;
  if(sGeometricModel.size() != 1)
  {
      std::cerr << "Unknown geometric model: " << sGeometricModel << std::endl;
      return EXIT_FAILURE;
  }
  const std::string sGeometricMode = std::string(1, std::tolower(sGeometricModel[0]));
  switch(sGeometricMode[0])
  {
    case 'f':
      eGeometricModelToCompute = FUNDAMENTAL_MATRIX;
    case 'e':
      eGeometricModelToCompute = ESSENTIAL_MATRIX;
    case 'h':
      eGeometricModelToCompute = HOMOGRAPHY_MATRIX;
    break;
    default:
      std::cerr << "Unknown geometric model: " << sGeometricMode << std::endl;
      return EXIT_FAILURE;
  }

  // -----------------------------
  // - Load SfM_Data Views & intrinsics data
  // a. Compute putative descriptor matches
  // b. Geometric filtering of putative matches
  // + Export some statistics
  // -----------------------------

  //---------------------------------------
  // Read SfM Scene (image view & intrinsics data)
  //---------------------------------------
  SfM_Data sfm_data;
  if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(VIEWS|INTRINSICS))) {
    std::cerr << std::endl
      << "The input SfM_Data file \""<< sSfM_Data_Filename << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }

  //---------------------------------------
  // Load SfM Scene regions
  //---------------------------------------
  // Init the regions_type from the image describer file (used for image regions extraction)
  using namespace openMVG::features;
  const std::string sImage_describer = stlplus::create_filespec(sMatchesDirectory, "image_describer", "json");
  std::unique_ptr<Regions> regions_type = Init_region_type_from_file(sImage_describer);
  if (!regions_type)
  {
    std::cerr << "Invalid: "
      << sImage_describer << " regions type file." << std::endl;
    return EXIT_FAILURE;
  }

  //---------------------------------------
  // a. Compute putative descriptor matches
  //    - Descriptor matching (according user method choice)
  //    - Keep correspondences only if NearestNeighbor ratio is ok
  //---------------------------------------

  // Load the corresponding view regions
  std::shared_ptr<Regions_Provider> regions_provider = std::make_shared<Regions_Provider>();
  if (!regions_provider->load(sfm_data, sMatchesDirectory, regions_type)) {
    std::cerr << std::endl << "Invalid regions." << std::endl;
    return EXIT_FAILURE;
  }

  PairWiseMatches map_PutativesMatches;

  // Build some alias from SfM_Data Views data:
  // - List views as a vector of filenames & image sizes
  std::vector<std::string> vec_fileNames;
  std::vector<std::pair<size_t, size_t> > vec_imagesSize;
  {
    vec_fileNames.reserve(sfm_data.GetViews().size());
    vec_imagesSize.reserve(sfm_data.GetViews().size());
    for (Views::const_iterator iter = sfm_data.GetViews().begin();
      iter != sfm_data.GetViews().end();
      ++iter)
    {
      const View * v = iter->second.get();
      vec_fileNames.push_back(stlplus::create_filespec(sfm_data.s_root_path,
          std::to_string(v->id_view)));
      vec_imagesSize.push_back( std::make_pair( v->ui_width, v->ui_height) );
    }
  }

  std::cout << std::endl << " - PUTATIVE MATCHES - " << std::endl;

  // If the matches already exists, reload them
  if (!bForce &&
     (!matchFilePerImage ?
          stlplus::is_file(stlplus::create_filespec(sMatchesDirectory, "putative.bin"))
        : stlplus::folder_wildcard(sMatchesDirectory, "*.putative.bin", false).size()))
  {
    Load(map_PutativesMatches, sfm_data.GetViewsKeys(), sMatchesDirectory, "putative");
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
    }

    // Allocate the right Matcher according the Matching requested method
    std::unique_ptr<Matcher> collectionMatcher;
    if (sNearestMatchingMethod == "AUTO")
    {
      if (regions_type->IsScalar())
      {
        std::cout << "Using FAST_CASCADE_HASHING_L2 matcher" << std::endl;
        collectionMatcher.reset(new Cascade_Hashing_Matcher_Regions_AllInMemory(fDistRatio));
      }
      else
      if (regions_type->IsBinary())
      {
        std::cout << "Using BRUTE_FORCE_HAMMING matcher" << std::endl;
        collectionMatcher.reset(new Matcher_Regions_AllInMemory(fDistRatio, BRUTE_FORCE_HAMMING));
      }
    }
    else
    if (sNearestMatchingMethod == "BRUTEFORCEL2")
    {
      std::cout << "Using BRUTE_FORCE_L2 matcher" << std::endl;
      collectionMatcher.reset(new Matcher_Regions_AllInMemory(fDistRatio, BRUTE_FORCE_L2));
    }
    else
    if (sNearestMatchingMethod == "BRUTEFORCEHAMMING")
    {
      std::cout << "Using BRUTE_FORCE_HAMMING matcher" << std::endl;
      collectionMatcher.reset(new Matcher_Regions_AllInMemory(fDistRatio, BRUTE_FORCE_HAMMING));
    }
    else
    if (sNearestMatchingMethod == "ANNL2")
    {
      std::cout << "Using ANN_L2 matcher" << std::endl;
      collectionMatcher.reset(new Matcher_Regions_AllInMemory(fDistRatio, ANN_L2));
    }
    else
    if (sNearestMatchingMethod == "CASCADEHASHINGL2")
    {
      std::cout << "Using CASCADE_HASHING_L2 matcher" << std::endl;
      collectionMatcher.reset(new Matcher_Regions_AllInMemory(fDistRatio, CASCADE_HASHING_L2));
    }
    else
    if (sNearestMatchingMethod == "FASTCASCADEHASHINGL2")
    {
      std::cout << "Using FAST_CASCADE_HASHING_L2 matcher" << std::endl;
      collectionMatcher.reset(new Cascade_Hashing_Matcher_Regions_AllInMemory(fDistRatio));
    }
    if (!collectionMatcher)
    {
      std::cerr << "Invalid Nearest Neighbor method: " << sNearestMatchingMethod << std::endl;
      return EXIT_FAILURE;
    }
    // Perform the matching
    system::Timer timer;
    {
      // From matching mode compute the pair list that have to be matched:
      Pair_Set pairs;
      switch (ePairmode)
      {
        case PAIR_EXHAUSTIVE: pairs = exhaustivePairs(sfm_data.GetViews(), rangeStart, rangeSize); break;
        case PAIR_CONTIGUOUS: pairs = contiguousWithOverlap(sfm_data.GetViews(), iMatchingVideoMode); break;
        case PAIR_FROM_FILE:
          std::cout << "Load pairList from file: " << sPredefinedPairList << std::endl;
          if(!loadPairs(sPredefinedPairList, pairs, orderPairs, rangeStart, rangeSize))
          {
              return EXIT_FAILURE;
          }
          break;
      }
      
      if( pairs.empty() )
      {
        std::cout << "No image pair to match." << std::endl;
        // If we only compute a selection of matches, we may have no match.
        return rangeSize ? EXIT_SUCCESS : EXIT_FAILURE;
      }
      std::cout << "There are " << sfm_data.GetViews().size() << " views and " << pairs.size() << " image pairs." << std::endl;

      // Photometric matching of putative pairs
      collectionMatcher->Match(sfm_data, regions_provider, pairs, map_PutativesMatches);
      
      if( map_PutativesMatches.empty() )
      {
        std::cout << "No putative matches." << std::endl;
        // If we only compute a selection of matches, we may have no match.
        return rangeSize ? EXIT_SUCCESS : EXIT_FAILURE;
      }
      std::cout << "There are " << map_PutativesMatches.size() << " putative matches." << std::endl;

      //---------------------------------------
      //-- Export putative matches
      //---------------------------------------
      Save(map_PutativesMatches, sMatchesDirectory, "putative", "bin", matchFilePerImage);
    }
    std::cout << "Task (Regions Matching) done in (s): " << timer.elapsed() << std::endl;
  }
  //-- export putative matches Adjacency matrix
  PairWiseMatchingToAdjacencyMatrixSVG(vec_fileNames.size(),
    map_PutativesMatches,
    stlplus::create_filespec(sMatchesDirectory, "PutativeAdjacencyMatrix", "svg"));
  //-- export view pair graph once putative graph matches have been computed
  {
    std::set<IndexT> set_ViewIds;
    std::transform(sfm_data.GetViews().begin(), sfm_data.GetViews().end(),
      std::inserter(set_ViewIds, set_ViewIds.begin()), stl::RetrieveKey());
    graph::indexedGraph putativeGraph(set_ViewIds, getPairs(map_PutativesMatches));
    graph::exportToGraphvizData(
      stlplus::create_filespec(sMatchesDirectory, "putative_matches.dot"),
      putativeGraph.g);
  }

  //---------------------------------------
  // b. Geometric filtering of putative matches
  //    - AContrario Estimation of the desired geometric model
  //    - Use an upper bound for the a contrario estimated threshold
  //---------------------------------------

  std::unique_ptr<ImageCollectionGeometricFilter> filter_ptr(
    new ImageCollectionGeometricFilter(&sfm_data, regions_provider));

  if (filter_ptr)
  {
    system::Timer timer;
    std::cout << std::endl << " - Geometric filtering - " << std::endl;

    PairWiseMatches map_GeometricMatches;
    switch (eGeometricModelToCompute)
    {
      case HOMOGRAPHY_MATRIX:
      {
        const bool bGeometric_only_guided_matching = true;
        filter_ptr->Robust_model_estimation(GeometricFilter_HMatrix_AC(4.0, imax_iteration),
          map_PutativesMatches, bGuided_matching,
          bGeometric_only_guided_matching ? -1.0 : 0.6);
        map_GeometricMatches = filter_ptr->Get_geometric_matches();
      }
      break;
      case FUNDAMENTAL_MATRIX:
      {
        filter_ptr->Robust_model_estimation(GeometricFilter_FMatrix_AC(4.0, imax_iteration),
          map_PutativesMatches, bGuided_matching);
        map_GeometricMatches = filter_ptr->Get_geometric_matches();
      }
      break;
      case ESSENTIAL_MATRIX:
      {
        filter_ptr->Robust_model_estimation(GeometricFilter_EMatrix_AC(4.0, imax_iteration),
          map_PutativesMatches, bGuided_matching);
        map_GeometricMatches = filter_ptr->Get_geometric_matches();

        //-- Perform an additional check to remove pairs with poor overlap
        std::vector<PairWiseMatches::key_type> vec_toRemove;
        for (PairWiseMatches::const_iterator iterMap = map_GeometricMatches.begin();
          iterMap != map_GeometricMatches.end(); ++iterMap)
        {
          const size_t putativePhotometricCount = map_PutativesMatches.find(iterMap->first)->second.size();
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
    }

    std::cout << "There are " << map_GeometricMatches.size() << " geometric matches." << std::endl;
    
    //---------------------------------------
    //-- Export geometric filtered matches
    //---------------------------------------
    Save(map_GeometricMatches, sMatchesDirectory, sGeometricMode, "bin", matchFilePerImage);

    std::cout << "Task done in (s): " << timer.elapsed() << std::endl;

    //-- export Adjacency matrix
    std::cout << "\n Export Adjacency Matrix of the pairwise's geometric matches"
      << std::endl;
    PairWiseMatchingToAdjacencyMatrixSVG(vec_fileNames.size(),
      map_GeometricMatches,
      stlplus::create_filespec(sMatchesDirectory, "GeometricAdjacencyMatrix", "svg"));

    //-- export view pair graph once geometric filter have been done
    {
      std::set<IndexT> set_ViewIds;
      std::transform(sfm_data.GetViews().begin(), sfm_data.GetViews().end(),
        std::inserter(set_ViewIds, set_ViewIds.begin()), stl::RetrieveKey());
      graph::indexedGraph putativeGraph(set_ViewIds, getPairs(map_GeometricMatches));
      graph::exportToGraphvizData(
        stlplus::create_filespec(sMatchesDirectory, "geometric_matches.dot"),
        putativeGraph.g);
    }
  }
  return EXIT_SUCCESS;
}
