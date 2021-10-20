// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2019 Pierre MOULON, Romuald Perrot.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/graph/graph.hpp"
#include "openMVG/graph/graph_stats.hpp"
#include "openMVG/matching/indMatch.hpp"
#include "openMVG/matching/indMatch_utils.hpp"
#include "openMVG/matching/pairwiseAdjacencyDisplay.hpp"
#include "openMVG/matching_image_collection/Cascade_Hashing_Matcher_Regions.hpp"
#include "openMVG/matching_image_collection/Matcher_Regions.hpp"
#include "openMVG/matching_image_collection/Pair_Builder.hpp"
#include "openMVG/sfm/pipelines/sfm_features_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_preemptive_regions_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_regions_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_regions_provider_cache.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/stl/stl.hpp"
#include "openMVG/system/timer.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>

using namespace openMVG;
using namespace openMVG::matching;
using namespace openMVG::sfm;
using namespace openMVG::matching_image_collection;
using namespace std;

/// Compute corresponding features between a series of views:
/// - Load view images description (regions: features & descriptors)
/// - Compute putative local feature matches (descriptors matching)
int main( int argc, char** argv )
{
  CmdLine cmd;

  std::string  sSfM_Data_Filename;
  std::string  sOutputMatchesFilename = "";
  float        fDistRatio             = 0.8f;
  std::string  sPredefinedPairList    = "";
  std::string  sNearestMatchingMethod = "AUTO";
  bool         bForce                 = false;
  unsigned int ui_max_cache_size      = 0;

  // Pre-emptive matching parameters
  unsigned int ui_preemptive_feature_count = 200;
  double preemptive_matching_percentage_threshold = 0.08;

  //required
  cmd.add( make_option( 'i', sSfM_Data_Filename, "input_file" ) );
  cmd.add( make_option( 'o', sOutputMatchesFilename, "out_dir" ) );
  cmd.add( make_option( 'p', sPredefinedPairList, "pair_list" ) );
  // Options
  cmd.add( make_option( 'r', fDistRatio, "ratio" ) );
  cmd.add( make_option( 'n', sNearestMatchingMethod, "nearest_matching_method" ) );
  cmd.add( make_option( 'f', bForce, "force" ) );
  cmd.add( make_option( 'c', ui_max_cache_size, "cache_size" ) );
  // Pre-emptive matching
  cmd.add( make_option( 'P', ui_preemptive_feature_count, "preemptive_feature_count") );


  try
  {
    if ( argc == 1 )
      throw std::string( "Invalid command line parameter." );
    cmd.process( argc, argv );
  }
  catch ( const std::string& s )
  {
    OPENMVG_LOG_INFO
      << "Usage: " << argv[ 0 ] << '\n'
      << "[-i|--input_file]   A SfM_Data file\n"
      << "[-o|--output_file]  Output file where computed matches are stored\n"
      << "[-p|--pair_list]    Pairs list file\n"
      << "\n[Optional]\n"
      << "[-f|--force] Force to recompute data]\n"
      << "[-r|--ratio] Distance ratio to discard non meaningful matches\n"
      << "   0.8: (default).\n"
      << "[-n|--nearest_matching_method]\n"
      << "  AUTO: auto choice from regions type,\n"
      << "  For Scalar based regions descriptor:\n"
      << "    BRUTEFORCEL2: L2 BruteForce matching,\n"
      << "    HNSWL2: L2 Approximate Matching with Hierarchical Navigable Small World graphs,\n"
      << "    HNSWL1: L1 Approximate Matching with Hierarchical Navigable Small World graphs\n"
      << "      tailored for quantized and histogram based descriptors (e.g uint8 RootSIFT)\n"
      << "    ANNL2: L2 Approximate Nearest Neighbor matching,\n"
      << "    CASCADEHASHINGL2: L2 Cascade Hashing matching.\n"
      << "    FASTCASCADEHASHINGL2: (default)\n"
      << "      L2 Cascade Hashing with precomputed hashed regions\n"
      << "     (faster than CASCADEHASHINGL2 but use more memory).\n"
      << "  For Binary based descriptor:\n"
      << "    BRUTEFORCEHAMMING: BruteForce Hamming matching,\n"
      << "    HNSWHAMMING: Hamming Approximate Matching with Hierarchical Navigable Small World graphs\n"
      << "[-c|--cache_size]\n"
      << "  Use a regions cache (only cache_size regions will be stored in memory)\n"
      << "  If not used, all regions will be load in memory."
      << "\n[Pre-emptive matching:]\n"
      << "[-P|--preemptive_feature_count] <NUMBER> Number of feature used for pre-emptive matching";

    OPENMVG_LOG_INFO << s;
    return EXIT_FAILURE;
  }

  OPENMVG_LOG_INFO << " You called : "
            << "\n"
            << argv[ 0 ] << "\n"
            << "--input_file " << sSfM_Data_Filename << "\n"
            << "--output_file " << sOutputMatchesFilename << "\n"
            << "--pair_list " << sPredefinedPairList << "\n"
            << "Optional parameters:"
            << "\n"
            << "--force " << bForce << "\n"
            << "--ratio " << fDistRatio << "\n"
            << "--nearest_matching_method " << sNearestMatchingMethod << "\n"
            << "--cache_size " << ((ui_max_cache_size == 0) ? "unlimited" : std::to_string(ui_max_cache_size)) << "\n"
            << "--preemptive_feature_used/count " << cmd.used('P') << " / " << ui_preemptive_feature_count;
  if (cmd.used('P'))
  {
    OPENMVG_LOG_INFO << "--preemptive_feature_count " << ui_preemptive_feature_count;
  }

  if ( sOutputMatchesFilename.empty() )
  {
    OPENMVG_LOG_ERROR << "No output file set.";
    return EXIT_FAILURE;
  }

  // -----------------------------
  // . Load SfM_Data Views & intrinsics data
  // . Compute putative descriptor matches
  // + Export some statistics
  // -----------------------------

  //---------------------------------------
  // Read SfM Scene (image view & intrinsics data)
  //---------------------------------------
  SfM_Data sfm_data;
  if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(VIEWS|INTRINSICS))) {
    OPENMVG_LOG_ERROR << "The input SfM_Data file \""<< sSfM_Data_Filename << "\" cannot be read.";
    return EXIT_FAILURE;
  }
  const std::string sMatchesDirectory = stlplus::folder_part( sOutputMatchesFilename );

  //---------------------------------------
  // Load SfM Scene regions
  //---------------------------------------
  // Init the regions_type from the image describer file (used for image regions extraction)
  using namespace openMVG::features;
  const std::string sImage_describer = stlplus::create_filespec(sMatchesDirectory, "image_describer", "json");
  std::unique_ptr<Regions> regions_type = Init_region_type_from_file(sImage_describer);
  if (!regions_type)
  {
    OPENMVG_LOG_ERROR << "Invalid: " << sImage_describer << " regions type file.";
    return EXIT_FAILURE;
  }

  //---------------------------------------
  // a. Compute putative descriptor matches
  //    - Descriptor matching (according user method choice)
  //    - Keep correspondences only if NearestNeighbor ratio is ok
  //---------------------------------------

  // Load the corresponding view regions
  std::shared_ptr<Regions_Provider> regions_provider;
  if (ui_max_cache_size == 0)
  {
    // Default regions provider (load & store all regions in memory)
    regions_provider = std::make_shared<Regions_Provider>();
  }
  else
  {
    // Cached regions provider (load & store regions on demand)
    regions_provider = std::make_shared<Regions_Provider_Cache>(ui_max_cache_size);
  }
  // If we use pre-emptive matching, we load less regions:
  if (ui_preemptive_feature_count > 0 && cmd.used('P'))
  {
    regions_provider = std::make_shared<Preemptive_Regions_Provider>(ui_preemptive_feature_count);
  }

  // Show the progress on the command line:
  system::LoggerProgress progress;

  if (!regions_provider->load(sfm_data, sMatchesDirectory, regions_type, &progress)) {
    OPENMVG_LOG_ERROR << "Cannot load view regions from: " << sMatchesDirectory << ".";
    return EXIT_FAILURE;
  }

  PairWiseMatches map_PutativeMatches;

  // Build some alias from SfM_Data Views data:
  // - List views as a vector of filenames & image sizes
  std::vector<std::string>               vec_fileNames;
  std::vector<std::pair<size_t, size_t>> vec_imagesSize;
  {
    vec_fileNames.reserve(sfm_data.GetViews().size());
    vec_imagesSize.reserve(sfm_data.GetViews().size());
    for (const auto view_it : sfm_data.GetViews())
    {
      const View * v = view_it.second.get();
      vec_fileNames.emplace_back(stlplus::create_filespec(sfm_data.s_root_path,
          v->s_Img_path));
      vec_imagesSize.emplace_back(v->ui_width, v->ui_height);
    }
  }

  OPENMVG_LOG_INFO << " - PUTATIVE MATCHES - ";
  // If the matches already exists, reload them
  if ( !bForce && ( stlplus::file_exists( sOutputMatchesFilename ) ) )
  {
    if ( !( Load( map_PutativeMatches, sOutputMatchesFilename ) ) )
    {
      OPENMVG_LOG_ERROR << "Cannot load input matches file";
      return EXIT_FAILURE;
    }
    OPENMVG_LOG_INFO
      << "\t PREVIOUS RESULTS LOADED;"
      << " #pair: " << map_PutativeMatches.size();
  }
  else // Compute the putative matches
  {
    // Allocate the right Matcher according the Matching requested method
    std::unique_ptr<Matcher> collectionMatcher;
    if ( sNearestMatchingMethod == "AUTO" )
    {
      if ( regions_type->IsScalar() )
      {
        OPENMVG_LOG_INFO << "Using FAST_CASCADE_HASHING_L2 matcher";
        collectionMatcher.reset(new Cascade_Hashing_Matcher_Regions(fDistRatio));
      }
      else
      if (regions_type->IsBinary())
      {
        OPENMVG_LOG_INFO << "Using HNSWHAMMING matcher";
        collectionMatcher.reset(new Matcher_Regions(fDistRatio, HNSW_HAMMING));
      }
    }
    else
    if (sNearestMatchingMethod == "BRUTEFORCEL2")
    {
      OPENMVG_LOG_INFO << "Using BRUTE_FORCE_L2 matcher";
      collectionMatcher.reset(new Matcher_Regions(fDistRatio, BRUTE_FORCE_L2));
    }
    else
    if (sNearestMatchingMethod == "BRUTEFORCEHAMMING")
    {
      OPENMVG_LOG_INFO << "Using BRUTE_FORCE_HAMMING matcher";
      collectionMatcher.reset(new Matcher_Regions(fDistRatio, BRUTE_FORCE_HAMMING));
    }
    else
    if (sNearestMatchingMethod == "HNSWL2")
    {
      OPENMVG_LOG_INFO << "Using HNSWL2 matcher";
      collectionMatcher.reset(new Matcher_Regions(fDistRatio, HNSW_L2));
    }
    if (sNearestMatchingMethod == "HNSWL1")
    {
      OPENMVG_LOG_INFO << "Using HNSWL1 matcher";
      collectionMatcher.reset(new Matcher_Regions(fDistRatio, HNSW_L1));
    }
    else
    if (sNearestMatchingMethod == "HNSWHAMMING")
    {
      OPENMVG_LOG_INFO << "Using HNSWHAMMING matcher";
      collectionMatcher.reset(new Matcher_Regions(fDistRatio, HNSW_HAMMING));
    }
    else
    if (sNearestMatchingMethod == "ANNL2")
    {
      OPENMVG_LOG_INFO << "Using ANN_L2 matcher";
      collectionMatcher.reset(new Matcher_Regions(fDistRatio, ANN_L2));
    }
    else
    if (sNearestMatchingMethod == "CASCADEHASHINGL2")
    {
      OPENMVG_LOG_INFO << "Using CASCADE_HASHING_L2 matcher";
      collectionMatcher.reset(new Matcher_Regions(fDistRatio, CASCADE_HASHING_L2));
    }
    else
    if (sNearestMatchingMethod == "FASTCASCADEHASHINGL2")
    {
      OPENMVG_LOG_INFO << "Using FAST_CASCADE_HASHING_L2 matcher";
      collectionMatcher.reset(new Cascade_Hashing_Matcher_Regions(fDistRatio));
    }
    if (!collectionMatcher)
    {
      OPENMVG_LOG_ERROR << "Invalid Nearest Neighbor method: " << sNearestMatchingMethod;
      return EXIT_FAILURE;
    }
    // Perform the matching
    system::Timer timer;
    {
      // From matching mode compute the pair list that have to be matched:
      Pair_Set pairs;
      if ( sPredefinedPairList.empty() )
      {
        OPENMVG_LOG_INFO << "No input pair file set. Use exhaustive match by default.";
        const size_t NImage = sfm_data.GetViews().size();
        pairs = exhaustivePairs( NImage );
      }
      else
      if ( !loadPairs( sfm_data.GetViews().size(), sPredefinedPairList, pairs ) )
      {
        OPENMVG_LOG_ERROR << "Failed to load pairs from file: \"" << sPredefinedPairList << "\"";
        return EXIT_FAILURE;
      }
      OPENMVG_LOG_INFO << "Running matching on #pairs: " << pairs.size();
      // Photometric matching of putative pairs
      collectionMatcher->Match( regions_provider, pairs, map_PutativeMatches, &progress );

      if (cmd.used('P')) // Preemptive filter
      {
        // Keep putative matches only if there is more than X matches
        PairWiseMatches map_filtered_matches;
        for (const auto & pairwisematches_it : map_PutativeMatches)
        {
          const size_t putative_match_count = pairwisematches_it.second.size();
          const int match_count_threshold =
            preemptive_matching_percentage_threshold * ui_preemptive_feature_count;
          // TODO: Add an option to keeping X Best pairs
          if (putative_match_count >= match_count_threshold)  {
            // the pair will be kept
            map_filtered_matches.insert(pairwisematches_it);
          }
        }
        map_PutativeMatches.clear();
        std::swap(map_filtered_matches, map_PutativeMatches);
      }

      //---------------------------------------
      //-- Export putative matches & pairs
      //---------------------------------------
      if ( !Save( map_PutativeMatches, std::string( sOutputMatchesFilename ) ) )
      {
        OPENMVG_LOG_ERROR
          << "Cannot save computed matches in: "
          << sOutputMatchesFilename;
        return EXIT_FAILURE;
      }
      // Save pairs
      const std::string sOutputPairFilename =
        stlplus::create_filespec( sMatchesDirectory, "preemptive_pairs", "txt" );
      if (!savePairs(
        sOutputPairFilename,
        getPairs(map_PutativeMatches)))
      {
        OPENMVG_LOG_ERROR
          << "Cannot save computed matches pairs in: "
          << sOutputPairFilename;
        return EXIT_FAILURE;
      }
    }
    OPENMVG_LOG_INFO << "Task (Regions Matching) done in (s): " << timer.elapsed();
  }

  OPENMVG_LOG_INFO << "#Putative pairs: " << map_PutativeMatches.size();

  // -- export Putative View Graph statistics
  graph::getGraphStatistics(sfm_data.GetViews().size(), getPairs(map_PutativeMatches));

  //-- export putative matches Adjacency matrix
  PairWiseMatchingToAdjacencyMatrixSVG( vec_fileNames.size(),
                                        map_PutativeMatches,
                                        stlplus::create_filespec( sMatchesDirectory, "PutativeAdjacencyMatrix", "svg" ) );
  //-- export view pair graph once putative graph matches has been computed
  {
    std::set<IndexT> set_ViewIds;
    std::transform( sfm_data.GetViews().begin(), sfm_data.GetViews().end(), std::inserter( set_ViewIds, set_ViewIds.begin() ), stl::RetrieveKey() );
    graph::indexedGraph putativeGraph( set_ViewIds, getPairs( map_PutativeMatches ) );
    graph::exportToGraphvizData(
        stlplus::create_filespec( sMatchesDirectory, "putative_matches" ),
        putativeGraph );
  }

  return EXIT_SUCCESS;
}
