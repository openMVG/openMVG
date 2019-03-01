// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2019 Pierre MOULON, Romuald PERROT

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/graph/graph.hpp"
#include "openMVG/features/akaze/image_describer_akaze.hpp"
#include "openMVG/features/descriptor.hpp"
#include "openMVG/features/feature.hpp"
#include "openMVG/matching/indMatch.hpp"
#include "openMVG/matching/indMatch_utils.hpp"
#include "openMVG/matching_image_collection/Matcher_Regions.hpp"
#include "openMVG/matching_image_collection/Cascade_Hashing_Matcher_Regions.hpp"
#include "openMVG/matching_image_collection/GeometricFilter.hpp"
#include "openMVG/sfm/pipelines/sfm_features_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_regions_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_regions_provider_cache.hpp"
#include "openMVG/matching_image_collection/F_ACRobust.hpp"
#include "openMVG/matching_image_collection/E_ACRobust.hpp"
#include "openMVG/matching_image_collection/E_ACRobust_Angular.hpp"
#include "openMVG/matching_image_collection/Eo_Robust.hpp"
#include "openMVG/matching_image_collection/H_ACRobust.hpp"
#include "openMVG/matching_image_collection/Pair_Builder.hpp"
#include "openMVG/matching/pairwiseAdjacencyDisplay.hpp"
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
using namespace openMVG::robust;
using namespace openMVG::sfm;
using namespace openMVG::matching_image_collection;
using namespace std;

enum EGeometricModel
{
  FUNDAMENTAL_MATRIX = 0,
  ESSENTIAL_MATRIX   = 1,
  HOMOGRAPHY_MATRIX  = 2,
  ESSENTIAL_MATRIX_ANGULAR = 3,
  ESSENTIAL_MATRIX_ORTHO = 4
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
  std::string sPutativePairsFilename;
  std::string sFilteredPairsFilename;
  std::string sGeometricModel = "f";
  bool bForce = false;
  bool bGuided_matching = false;
  int imax_iteration = 2048;
  unsigned int ui_max_cache_size = 0;

  //required
  cmd.add( make_option('i', sSfM_Data_Filename, "input_file") );
  cmd.add( make_option('o', sFilteredPairsFilename, "output_file") );
  cmd.add( make_option('p', sPutativePairsFilename, "pairs") );
  // Options
  cmd.add( make_option('g', sGeometricModel, "geometric_model") );
  cmd.add( make_option('f', bForce, "force") );
  cmd.add( make_option('m', bGuided_matching, "guided_matching") );
  cmd.add( make_option('I', imax_iteration, "max_iteration") );
  cmd.add( make_option('c', ui_max_cache_size, "cache_size") );


  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch (const std::string& s) {
      std::cerr << "Usage: " << argv[0] << '\n'
      << "[-i|--input_file] a SfM_Data file\n"
      << "[-p|--pairs] initial pairs filename to use\n"
      << "[-o|--output_file] output filtered pairs filename\n"
      << "\n[Optional]\n"
      << "[-f|--force] Force to recompute data]\n"
      << "[-g|--geometric_model]\n"
      << "  (pairwise correspondences filtering thanks to robust model estimation):\n"
      << "   f: (default) fundamental matrix,\n"
      << "   e: essential matrix,\n"
      << "   h: homography matrix.\n"
      << "   a: essential matrix with an angular parametrization,\n"
      << "   o: orthographic essential matrix.\n"
      << "[-m|--guided_matching]\n"
      << "  use the found model to improve the pairwise correspondences.\n"
      << "[-c|--cache_size]\n"
      << "  Use a regions cache (only cache_size regions will be stored in memory)\n"
      << "  If not used, all regions will be load in memory."
      << std::endl;

      std::cerr << s << std::endl;
      return EXIT_FAILURE;
  }

  std::cout << " You called : " << "\n"
            << argv[0] << "\n"
            << "--input_file: " << sSfM_Data_Filename << "\n"
            << "--pairs: " << sPutativePairsFilename << "\n"
            << "--output_file: " << sFilteredPairsFilename << "\n"
            << "Optional parameters:" << "\n"
            << "--force " << bForce << "\n"
            << "--geometric_model " << sGeometricModel << "\n"
            << "--guided_matching " << bGuided_matching << "\n"
            << "--cache_size " << ((ui_max_cache_size == 0) ? "unlimited" : std::to_string(ui_max_cache_size)) << std::endl;



  if (sFilteredPairsFilename.empty() )  
  {
    std::cerr << "\nIt is an invalid output file" << std::endl;
    return EXIT_FAILURE;
  }
  if( sSfM_Data_Filename.empty() )
  {
    std::cerr << "\nIt is an invalid SfM file" << std::endl; 
    return EXIT_FAILURE;
  }
  if( sPutativePairsFilename.empty() )
  {
    std::cerr << "\nIt is an invalid putative pairs file" << std::endl; 
    return EXIT_FAILURE;
  }

  const std::string sMatchesDirectory = stlplus::folder_part( sSfM_Data_Filename ) + "/../matches/" ;


  EGeometricModel eGeometricModelToCompute = FUNDAMENTAL_MATRIX;
  switch (sGeometricModel[0])
  {
    case 'f': case 'F':
      eGeometricModelToCompute = FUNDAMENTAL_MATRIX;
    break;
    case 'e': case 'E':
      eGeometricModelToCompute = ESSENTIAL_MATRIX;
    break;
    case 'h': case 'H':
      eGeometricModelToCompute = HOMOGRAPHY_MATRIX;
    break;
    case 'a': case 'A':
      eGeometricModelToCompute = ESSENTIAL_MATRIX_ANGULAR;
    break;
    case 'o': case 'O':
      eGeometricModelToCompute = ESSENTIAL_MATRIX_ORTHO;
    break;
    default:
      std::cerr << "Unknown geometric model" << std::endl;
      return EXIT_FAILURE;
  }

  // -----------------------------
  // - Load SfM_Data Views & intrinsics data
  // a. Load putative descriptor matches
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
  // Consider that the image_describer.json is inside the matches directory (which is bellow the sfm_data.bin)
  const std::string sImage_describer = stlplus::create_filespec( sMatchesDirectory , "image_describer.json" ) ;
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

  // Show the progress on the command line:
  C_Progress_display progress;

  if (!regions_provider->load(sfm_data, sMatchesDirectory, regions_type, &progress)) 
  {
    std::cerr << std::endl << "Invalid regions." << std::endl;
    return EXIT_FAILURE;
  }

  PairWiseMatches map_PutativesMatches;
  //---------------------------------------
  // A. Load initial matches 
  //---------------------------------------
  if( ! Load( map_PutativesMatches, sPutativePairsFilename) ) 
  {
    std::cerr << "Filed to load the initial pairs" ;
    return EXIT_FAILURE;
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
    const double d_distance_ratio = 0.6;

    PairWiseMatches map_GeometricMatches;
    switch (eGeometricModelToCompute)
    {
      case HOMOGRAPHY_MATRIX:
      {
        const bool bGeometric_only_guided_matching = true;
        filter_ptr->Robust_model_estimation(
          GeometricFilter_HMatrix_AC(4.0, imax_iteration),
          map_PutativesMatches, bGuided_matching,
          bGeometric_only_guided_matching ? -1.0 : d_distance_ratio, &progress);
        map_GeometricMatches = filter_ptr->Get_geometric_matches();
      }
      break;
      case FUNDAMENTAL_MATRIX:
      {
        filter_ptr->Robust_model_estimation(
          GeometricFilter_FMatrix_AC(4.0, imax_iteration),
          map_PutativesMatches, bGuided_matching, d_distance_ratio, &progress);
        map_GeometricMatches = filter_ptr->Get_geometric_matches();
      }
      break;
      case ESSENTIAL_MATRIX:
      {
        filter_ptr->Robust_model_estimation(
          GeometricFilter_EMatrix_AC(4.0, imax_iteration),
          map_PutativesMatches, bGuided_matching, d_distance_ratio, &progress);
        map_GeometricMatches = filter_ptr->Get_geometric_matches();

        //-- Perform an additional check to remove pairs with poor overlap
        std::vector<PairWiseMatches::key_type> vec_toRemove;
        for (const auto & pairwisematches_it : map_GeometricMatches)
        {
          const size_t putativePhotometricCount = map_PutativesMatches.find(pairwisematches_it.first)->second.size();
          const size_t putativeGeometricCount = pairwisematches_it.second.size();
          const float ratio = putativeGeometricCount / static_cast<float>(putativePhotometricCount);
          if (putativeGeometricCount < 50 || ratio < .3f)  {
            // the pair will be removed
            vec_toRemove.push_back(pairwisematches_it.first);
          }
        }
        //-- remove discarded pairs
        for (const auto & pair_to_remove_it : vec_toRemove)
        {
          map_GeometricMatches.erase(pair_to_remove_it);
        }
      }
      break;
      case ESSENTIAL_MATRIX_ANGULAR:
      {
        filter_ptr->Robust_model_estimation(
          GeometricFilter_ESphericalMatrix_AC_Angular(4.0, imax_iteration),
          map_PutativesMatches, bGuided_matching);
        map_GeometricMatches = filter_ptr->Get_geometric_matches();
      }
      break;
      case ESSENTIAL_MATRIX_ORTHO:
      {
        filter_ptr->Robust_model_estimation(
          GeometricFilter_EOMatrix_RA(2.0, imax_iteration),
          map_PutativesMatches, bGuided_matching, d_distance_ratio, &progress);
        map_GeometricMatches = filter_ptr->Get_geometric_matches();
      }
      break;
    }

    //---------------------------------------
    //-- Export geometric filtered matches
    //---------------------------------------
    if (!Save(map_GeometricMatches, sFilteredPairsFilename) )
    {
      std::cerr
          << "Cannot save filtered matches in: " 
          << sFilteredPairsFilename ;
      return EXIT_FAILURE;
    }

    std::cout << "Task done in (s): " << timer.elapsed() << std::endl;

    //-- export Adjacency matrix
    std::cout << "\n Export Adjacency Matrix of the pairwise's geometric matches"
      << std::endl;
    PairWiseMatchingToAdjacencyMatrixSVG(sfm_data.GetViews().size(),
      map_GeometricMatches,
      stlplus::create_filespec(sMatchesDirectory, "GeometricAdjacencyMatrix", "svg"));

    //-- export view pair graph once geometric filter have been done
    {
      std::set<IndexT> set_ViewIds;
      std::transform(sfm_data.GetViews().begin(), sfm_data.GetViews().end(),
        std::inserter(set_ViewIds, set_ViewIds.begin()), stl::RetrieveKey());
      graph::indexedGraph putativeGraph(set_ViewIds, getPairs(map_GeometricMatches));
      graph::exportToGraphvizData(
        stlplus::create_filespec(sMatchesDirectory, "geometric_matches"),
        putativeGraph);
    }
  }
  return EXIT_SUCCESS;
}
