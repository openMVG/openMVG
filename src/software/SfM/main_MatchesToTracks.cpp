// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2018 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/features/feature.hpp"
#include "openMVG/image/image_io.hpp"
#include "openMVG/matching/indMatch.hpp"
#include "openMVG/matching/indMatch_utils.hpp"
#include "openMVG/matching/svg_matches.hpp"
#include "openMVG/sfm/pipelines/sfm_features_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_matches_provider.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/tracks/tracks.hpp"

#include "software/SfM/SfMIOHelper.hpp"
#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "third_party/progress/progress_display.hpp"

using namespace openMVG;
using namespace openMVG::matching;
using namespace openMVG::sfm;
using namespace openMVG::tracks;

int main(int argc, char ** argv)
{
  CmdLine cmd;

  std::string sSfM_Data_Filename;
  std::string sMatchesDir;
  std::string sMatchFile;
  std::string sOutDir = "";

  cmd.add( make_option('i', sSfM_Data_Filename, "input_file") );
  cmd.add( make_option('d', sMatchesDir, "matchdir") );
  cmd.add( make_option('m', sMatchFile, "matchfile") );
  cmd.add( make_option('o', sOutDir, "outdir") );

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch (const std::string& s) {
      std::cerr << "Export pairwise tracks.\nUsage: " << argv[0] << "\n"
      << "[-i|--input_file file] path to a SfM_Data scene\n"
      << "[-d|--matchdir path]\n"
      << "[-m|--sMatchFile filename]\n"
      << "[-o|--outdir path]\n"
      << std::endl;

      std::cerr << s << std::endl;
      return EXIT_FAILURE;
  }

  if (sOutDir.empty())  {
    std::cerr << "\nIt is an invalid output directory" << std::endl;
    return EXIT_FAILURE;
  }


  //---------------------------------------
  // Read SfM Scene (image view names)
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
  const std::string sImage_describer = stlplus::create_filespec(sMatchesDir, "image_describer", "json");
  std::unique_ptr<Regions> regions_type = Init_region_type_from_file(sImage_describer);
  if (!regions_type)
  {
    std::cerr << "Invalid: " << sImage_describer << " regions type file." << std::endl;
    return EXIT_FAILURE;
  }

  // Read the features
  std::shared_ptr<Features_Provider> feats_provider = std::make_shared<Features_Provider>();
  if (!feats_provider->load(sfm_data, sMatchesDir, regions_type)) {
    std::cerr << "Invalid features." << std::endl;
    return EXIT_FAILURE;
  }

  // Read the matches
  std::shared_ptr<Matches_Provider> matches_provider = std::make_shared<Matches_Provider>();
  if (!matches_provider->load(sfm_data, sMatchFile)) {
    std::cerr << "Invalid match file." << std::endl;
    return EXIT_FAILURE;
  }

  //---------------------------------------
  // Compute tracks from matches
  //---------------------------------------
  tracks::STLMAPTracks map_tracks;
  {
    const openMVG::matching::PairWiseMatches & map_Matches = matches_provider->pairWise_matches_;
    tracks::TracksBuilder tracksBuilder;
    tracksBuilder.Build(map_Matches);
    tracksBuilder.Filter();
    tracksBuilder.ExportToSTL(map_tracks);
  }

  // Init the putative landmarks
  {
    // For every track add the obervations:
    // - views and feature positions that see this landmark
    for ( const auto & iterT : map_tracks )
    {
      Observations obs;
      for (const auto & track_ids : iterT.second) // {ViewId, FeatureId}
      {
        const auto & view_id = track_ids.first;
        const auto & feat_id = track_ids.second;
        const Vec2 x = feats_provider->feats_per_view[view_id][feat_id].coords().template cast<double>();
        obs.insert({view_id, Observation(x, feat_id)});
      }
      sfm_data.structure[iterT.first].obs = std::move(obs);
    }
  }

  if (!Save(sfm_data, stlplus::create_filespec(sOutDir, "sfm_data.bin"),ESfM_Data(ALL)))
  {
    std::cerr << "Unable to write the output sfm_data.bin file" << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
