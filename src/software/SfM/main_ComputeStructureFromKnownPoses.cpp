// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/cameras/Camera_Common.hpp"
#include "openMVG/features/feature.hpp"
#include "openMVG/features/svg_features.hpp"
#include "openMVG/geometry/frustum.hpp"
#include "openMVG/matching/indMatch.hpp"
#include "openMVG/matching/indMatch_utils.hpp"
#include "openMVG/matching_image_collection/Pair_Builder.hpp"
#include "openMVG/sfm/pipelines/structure_from_known_poses/structure_estimator.hpp"
#include "openMVG/sfm/pipelines/sfm_regions_provider_cache.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_BA.hpp"
#include "openMVG/sfm/sfm_data_BA_ceres.hpp"
#include "openMVG/sfm/sfm_data_filters.hpp"
#include "openMVG/sfm/sfm_data_filters_frustum.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_report.hpp"
#include "openMVG/sfm/sfm_data_triangulation.hpp"
#include "openMVG/tracks/tracks.hpp"
#include "openMVG/system/timer.hpp"
#include "openMVG/types.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/progress/progress_display.hpp"

#include <iostream>
#include <memory>
#include <string>

using namespace openMVG;
using namespace openMVG::matching;
using namespace openMVG::sfm;

/// Build a list of pair from the camera frusta intersections
Pair_Set BuildPairsFromFrustumsIntersections(
  const SfM_Data & sfm_data,
  const double z_near = -1., // default near plane
  const double z_far = -1.)  // default far plane
{
  const Frustum_Filter frustum_filter(sfm_data, z_near, z_far);
  return frustum_filter.getFrustumIntersectionPairs();
}

/// Compute the structure of a scene according existing camera poses.
int main(int argc, char **argv)
{
  using namespace std;
  std::cout << "Compute Structure from the provided poses" << std::endl;

  CmdLine cmd;

  std::string sSfM_Data_Filename;
  std::string sMatchesDir;
  std::string sMatchFile;
  std::string sPairFile;
  std::string sOutFile = "";
  double dMax_reprojection_error = 4.0;
  unsigned int ui_max_cache_size = 0;

  cmd.add( make_option('i', sSfM_Data_Filename, "input_file") );
  cmd.add( make_option('m', sMatchesDir, "match_dir") );
  cmd.add( make_option('f', sMatchFile, "match_file") );
  cmd.add( make_option('p', sPairFile, "pair_file") );
  cmd.add( make_option('o', sOutFile, "output_file") );
  cmd.add( make_switch('b', "bundle_adjustment"));
  cmd.add( make_option('r', dMax_reprojection_error, "residual_threshold"));
  cmd.add( make_option('c', ui_max_cache_size, "cache_size") );
  cmd.add( make_switch('d', "direct_triangulation"));

  try {
    if (argc == 1) throw std::string("Invalid command line parameter.");
    cmd.process(argc, argv);
  } catch (const std::string& s) {
    std::cerr << "Usage: " << argv[0] << '\n'
    << "[-i|--input_file] path to a SfM_Data scene\n"
    << "[-m|--match_dir] path to the features and descriptors that "
    <<    "corresponds to the provided SfM_Data scene\n"
    << "[-o|--output_file] file where the output data will be stored "
    <<    "(i.e. path/sfm_data_structure.bin)\n"

    << "\n[Triangulation mode]:\n"
    << " [No Provided Matches -> Triangulation of guided epipolar geometry matches (default mode)]\n"
    << "\t[-p|--pair_file] path to a pairs file (only those pairs will be considered to compute the structure)\n"
    << "\t[-f|--match_file] path to a matches file (loaded pair indexes will be used)\n"

    << " [Provided Matches -> Robust triangulation of the match file (activated by -d)]\n"
    << "\t[-d|--direct_triangulation] Robustly triangulate the tracks computed from the file given by [--f|--match_file]\n"
    << "\t[-f|--match_file] path to a matches file (loaded pair indexes will be used)\n"

    << "\n[Optional]\n"
    << "[-b|--bundle_adjustment] (switch) perform a bundle adjustment on the scene (OFF by default)\n"
    << "[-r|--residual_threshold] maximal pixels reprojection error that will be considered for triangulations (4.0 by default)\n"
        << "[-c|--cache_size]\n"
    << "  Use a regions cache (only cache_size regions will be stored in memory)\n"
    << "  If not used, all regions will be load in memory.\n"

    << std::endl;

    std::cerr << s << std::endl;
    return EXIT_FAILURE;
  }

  // Load input SfM_Data scene
  SfM_Data sfm_data;
  if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(VIEWS|INTRINSICS|EXTRINSICS))) {
    std::cerr << std::endl
      << "The input SfM_Data file \""<< sSfM_Data_Filename << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }

  // Init the regions_type from the image describer file (used for image regions extraction)
  using namespace openMVG::features;
  const std::string sImage_describer = stlplus::create_filespec(sMatchesDir, "image_describer", "json");
  std::unique_ptr<Regions> regions_type = Init_region_type_from_file(sImage_describer);
  if (!regions_type)
  {
    std::cerr << "Invalid: "
      << sImage_describer << " regions type file." << std::endl;
    return EXIT_FAILURE;
  }

  // Prepare the Regions provider
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

  if (!regions_provider->load(sfm_data, sMatchesDir, regions_type, &progress)) {
    std::cerr << std::endl
      << "Invalid regions." << std::endl;
    return EXIT_FAILURE;
  }

  std::cout
    << "Loaded a sfm_data scene with:\n"
    << " #views: " << sfm_data.GetViews().size() << "\n"
    << " #poses: " << sfm_data.GetPoses().size() << "\n"
    << " #intrinsics: " << sfm_data.GetIntrinsics().size() <<  "\n"
    << " #tracks: " << sfm_data.GetLandmarks().size()
    << std::endl;

  const bool bDirect_triangulation = cmd.used('d');

  if (bDirect_triangulation)
  {
    // Check that a match file have been provided
    if (sMatchFile.empty() || !sPairFile.empty())
    {
      std::cerr << "You must provide a match file thanks to the [--f|--match_file] option" << std::endl;
      return EXIT_FAILURE;
    }
    std::cout
      << "\n======================================\n"
      << "Robust triangulation of the match file\n"
      << "======================================"  << std::endl;
    PairWiseMatches matches;
    if (!Load(matches, sMatchFile))
    {
      std::cerr<< "Unable to read the matches file." << std::endl;
      return EXIT_FAILURE;
    }
    // Compute the tracks from the pairwise estimation
    // Compute tracks from matches
    const int min_track_length = 2;
    openMVG::tracks::STLMAPTracks tracks;
    {
      // List of features matches for each couple of images
      std::cout << "\n" << "Building tracks..." << std::endl;
      tracks::TracksBuilder tracks_builder;
      tracks_builder.Build(matches);
      std::cout << "Filtering tracks..." << std::endl;
      tracks_builder.Filter(min_track_length);
      //-- Build tracks with STL compliant type :
      tracks_builder.ExportToSTL(tracks);

      // Display some statistics about the computed tracks
      {
        std::ostringstream track_stream;
        //-- Display stats :
        //    - number of images
        //    - number of tracks
        std::set<uint32_t> images_id;
        tracks::TracksUtilsMap::ImageIdInTracks(tracks, images_id);
        track_stream
          << "------------------" << "\n"
          << "-- Tracks Stats --" << "\n"
          << " Tracks number: " << tracks_builder.NbTracks() << "\n"
          << " Images Id: " << "\n";
        std::copy(images_id.begin(), images_id.end(),
          std::ostream_iterator<uint32_t>(track_stream, ", "));
        track_stream << "\n------------------" << "\n";

        std::map<uint32_t, uint32_t> track_length_histogram;
        tracks::TracksUtilsMap::TracksLength(tracks, track_length_histogram);
        track_stream << "TrackLength, Count" << "\n";
        for (const auto & it : track_length_histogram)  {
          track_stream << "\t" << it.first << "\t" << it.second << "\n";
        }
        track_stream << "\n";
        std::cout << track_stream.str();
      }
    }

    std::cout
      << "====================================\n"
      << "Robust triangulation of the tracks\n"
      << " - tracks computed from a match file\n"
      << "====================================" << std::endl;

    // Fill sfm_data with the computed tracks (no 3D yet)
    Landmarks & structure = sfm_data.structure;
    IndexT idx(0);
    for (const auto & tracks_it : tracks)
    {
      structure[idx] = {};
      Observations & obs = structure.at(idx).obs;
      for (const auto & track_it : tracks_it.second)
      {
        const auto imaIndex = track_it.first;
        const auto featIndex = track_it.second;
        const Vec2 & pt = regions_provider->get(imaIndex)->GetRegionPosition(featIndex);
        obs[imaIndex] = {pt, featIndex};
      }
      ++idx;
    }

    // Compute 3D position of the landmark of the structure by robust triangulation of the observations
    {
      const double max_reprojection_error = 4.0; // pixels reprojection error
      bool console_verbose = true;
      SfM_Data_Structure_Computation_Robust structure_estimator(
        max_reprojection_error,
        min_track_length,
        min_track_length,
        console_verbose);
      structure_estimator.triangulate(sfm_data);
    }
  }
  else
  {
    std::cout
      << "=============================================================\n"
      << "Robust triangulation of the tracks\n"
      << " - Triangulation of guided epipolar geometry matches\n"
      << "============================================================="
      << std::endl;
    //--
    //- Pair selection method:
    //  - geometry guided -> camera frustum intersection,
    //  - putative matches guided (photometric matches)
    //     (keep pairs that have valid Intrinsic & Pose ids).
    //--
    Pair_Set pairs;
    if (sMatchFile.empty() && sPairFile.empty())
    {
      // no provided pair, use camera frustum intersection
      pairs = BuildPairsFromFrustumsIntersections(sfm_data);
    }
    else
    {
      if (!sPairFile.empty() && sMatchFile.empty())
      {
        if (!loadPairs(sfm_data.GetViews().size(), sPairFile, pairs))
        {
          std::cerr << "Unable to read the pair file." << std::endl;
          return EXIT_FAILURE;
        }
      }
      else if (!sMatchFile.empty() && sPairFile.empty())
      {
        PairWiseMatches matches;
        if (!Load(matches, sMatchFile))
        {
          std::cerr<< "Unable to read the matches file." << std::endl;
          return EXIT_FAILURE;
        }
        pairs = getPairs(matches);
        // Keep only Pairs that belong to valid view indexes.
        const std::set<IndexT> valid_viewIdx = Get_Valid_Views(sfm_data);
        pairs = Pair_filter(pairs, valid_viewIdx);
      }
      else
      {
        std::cerr << "Cannot use --match_file and --pair_file at the same time" << std::endl;
      }
    }

    openMVG::system::Timer timer;

    //------------------------------------------
    // Compute Structure from known camera poses
    //------------------------------------------
    SfM_Data_Structure_Estimation_From_Known_Poses structure_estimator(dMax_reprojection_error);
    structure_estimator.run(sfm_data, pairs, regions_provider);
    std::cout << "\nStructure estimation took (s): " << timer.elapsed() << "." << std::endl;

  }
  regions_provider.reset(); // Regions are not longer needed.
  RemoveOutliers_AngleError(sfm_data, 2.0);

  std::cout
    << "\n#landmark found: " << sfm_data.GetLandmarks().size() << std::endl;

  std::cout << "...Generating SfM_Report.html" << std::endl;
  Generate_SfM_Report(sfm_data,
    stlplus::create_filespec(stlplus::folder_part(sOutFile), "SfMStructureFromKnownPoses_Report.html"));

  if (cmd.used('b'))
  {
    // Check that poses & intrinsic cover some measures (after outlier removal)
    const IndexT minPointPerPose = 12; // 6 min
    const IndexT minTrackLength = 3; // 2 min
    if (eraseUnstablePosesAndObservations(sfm_data, minPointPerPose, minTrackLength))
    {
      KeepLargestViewCCTracks(sfm_data);
      eraseUnstablePosesAndObservations(sfm_data, minPointPerPose, minTrackLength);

      const size_t pointcount_cleaning = sfm_data.structure.size();
      std::cout << "Point_cloud cleaning:\n"
        << "\t #3DPoints: " << pointcount_cleaning << "\n";
    }

    std::cout << "Bundle adjustment..." << std::endl;
    Bundle_Adjustment_Ceres bundle_adjustment_obj;
    bundle_adjustment_obj.Adjust
      (
        sfm_data,
        Optimize_Options(
          cameras::Intrinsic_Parameter_Type::ADJUST_ALL,
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL)
      );
  }

  std::cout
    << "Found a sfm_data scene with:\n"
    << " #views: " << sfm_data.GetViews().size() << "\n"
    << " #poses: " << sfm_data.GetPoses().size() << "\n"
    << " #intrinsics: " << sfm_data.GetIntrinsics().size() <<  "\n"
    << " #tracks: " << sfm_data.GetLandmarks().size()
    << std::endl;

  if (stlplus::extension_part(sOutFile) != "ply") {
    Save(sfm_data,
      stlplus::create_filespec(
        stlplus::folder_part(sOutFile),
        stlplus::basename_part(sOutFile), "ply"),
      ESfM_Data(ALL));
  }

  if (Save(sfm_data, sOutFile, ESfM_Data(ALL)))
  {
    return EXIT_SUCCESS;
  }
  return EXIT_FAILURE;
}
