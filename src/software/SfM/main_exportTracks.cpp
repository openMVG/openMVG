// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013, 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/matching/indMatch.hpp"
#include "openMVG/matching/indMatch_utils.hpp"
#include "openMVG/matching/svg_matches.hpp"
#include "openMVG/image/image_io.hpp"
#include "openMVG/features/feature.hpp"
#include "openMVG/tracks/tracks.hpp"
#include "openMVG/sfm/pipelines/sfm_features_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_matches_provider.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"

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
    std::cerr << "Invalid: "
      << sImage_describer << " regions type file." << std::endl;
    return EXIT_FAILURE;
  }

  // Read the features
  std::shared_ptr<Features_Provider> feats_provider = std::make_shared<Features_Provider>();
  if (!feats_provider->load(sfm_data, sMatchesDir, regions_type)) {
    std::cerr << std::endl
      << "Invalid features." << std::endl;
    return EXIT_FAILURE;
  }
  // Read the matches
  std::shared_ptr<Matches_Provider> matches_provider = std::make_shared<Matches_Provider>();
  if (!matches_provider->load(sfm_data, sMatchFile)) {
    std::cerr << "\nInvalid matches file." << std::endl;
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
  openMVG::tracks::SharedTrackVisibilityHelper track_visibility_helper(map_tracks);

  // ------------
  // For each pair, export the matches
  // ------------
  const uint32_t viewCount(sfm_data.GetViews().size());

  stlplus::folder_create(sOutDir);
  std::cout << "\n viewCount: " << viewCount << std::endl;
  std::cout << "\n Export pairwise tracks" << std::endl;
  C_Progress_display my_progress_bar( (viewCount*(viewCount-1)) / 2.0 );

  for (uint32_t I = 0; I < viewCount; ++I)
  {
    for (uint32_t J = I+1; J < viewCount; ++J, ++my_progress_bar)
    {
      const View
        *view_I = sfm_data.GetViews().at(I).get(),
        *view_J = sfm_data.GetViews().at(J).get();

      const std::string
        sView_I = stlplus::create_filespec(sfm_data.s_root_path, view_I->s_Img_path),
        sView_J = stlplus::create_filespec(sfm_data.s_root_path, view_J->s_Img_path);

      // Get common tracks between view I and J
      tracks::STLMAPTracks map_tracksCommon;
      track_visibility_helper.GetTracksInImages({I,J}, map_tracksCommon);

      if (!map_tracksCommon.empty())
      {
        // Build corresponding indexes from the two view tracks
        matching::IndMatches matches;
        matches.reserve(map_tracksCommon.size());
        for (const auto & tracks_it : map_tracksCommon)
        {
          tracks::submapTrack::const_iterator iter = tracks_it.second.begin();
          const IndexT i = iter->second; ++iter;
          const IndexT j = iter->second;
          matches.emplace_back(i, j);
        }

        // Draw corresponding features
        const bool bVertical = false;
        std::ostringstream os;
        os << stlplus::folder_append_separator(sOutDir)
           << I << "_" << J
           << "_" << map_tracksCommon.size() << "_.svg";
        Matches2SVG
        (
          sView_I,
          {view_I->ui_width, view_I->ui_height},
          feats_provider->getFeatures(view_I->id_view),
          sView_J,
          {view_J->ui_width, view_J->ui_height},
          feats_provider->getFeatures(view_J->id_view),
          matches,
          os.str(),
          bVertical
        );
      }
    }
  }
  return EXIT_SUCCESS;
}
