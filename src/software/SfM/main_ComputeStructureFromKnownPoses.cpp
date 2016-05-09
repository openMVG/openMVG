
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/sfm.hpp"
#include "openMVG/system/timer.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

using namespace openMVG;
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
  std::string sOutFile = "";
  double dMax_reprojection_error = 4.0;

  cmd.add( make_option('i', sSfM_Data_Filename, "input_file") );
  cmd.add( make_option('m', sMatchesDir, "match_dir") );
  cmd.add( make_option('f', sMatchFile, "match_file") );
  cmd.add( make_option('o', sOutFile, "output_file") );
  cmd.add( make_switch('b', "bundle_adjustment"));
  cmd.add( make_option('r', dMax_reprojection_error, "residual_threshold"));

  try {
    if (argc == 1) throw std::string("Invalid command line parameter.");
    cmd.process(argc, argv);
  } catch(const std::string& s) {
    std::cerr << "Usage: " << argv[0] << '\n'
    << "[-i|--input_file] path to a SfM_Data scene\n"
    << "[-m|--match_dir] path to the features and descriptors that "
    <<    "corresponds to the provided SfM_Data scene\n"
    << "[-o|--output_file] file where the output data will be stored "
    <<    "(i.e. path/sfm_data_structure.bin)\n"
    << "\n[Optional]\n"
    << "[-f|--match_file] path to a matches file (loaded pair indexes will be used)\n"
    << "[-b|--bundle_adjustment] (switch) perform a bundle adjustment on the scene (OFF by default)\n"
    << "[-r|--residual_threshold] maximal pixels reprojection error that will be considered for triangulations (4.0 by default)\n"
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
  std::shared_ptr<Regions_Provider> regions_provider = std::make_shared<Regions_Provider>();
  if (!regions_provider->load(sfm_data, sMatchesDir, regions_type)) {
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

  //--
  //- Pair selection method:
  //  - geometry guided -> camera frustum intersection,
  //  - putative matches guided (photometric matches)
  //     (keep pairs that have valid Intrinsic & Pose ids).
  //--
  Pair_Set pairs;
  if (sMatchFile.empty())
  {
    // no provided pair, use camera frustum intersection
    pairs = BuildPairsFromFrustumsIntersections(sfm_data);
  }
  else
  {
    PairWiseMatches matches;
    if (!matching::Load(matches, sMatchFile)) {
      std::cerr<< "Unable to read the matches file." << std::endl;
      return EXIT_FAILURE;
    }
    pairs = getPairs(matches);
    // Keep only Pairs that belong to valid view indexes.
    const std::set<IndexT> valid_viewIdx = Get_Valid_Views(sfm_data);
    pairs = Pair_filter(pairs, valid_viewIdx);
  }

  openMVG::system::Timer timer;

  //------------------------------------------
  // Compute Structure from known camera poses
  //------------------------------------------
  SfM_Data_Structure_Estimation_From_Known_Poses structure_estimator(dMax_reprojection_error);
  structure_estimator.run(sfm_data, pairs, regions_provider);
  regions_provider.reset(); // Regions are not longer needed.
  RemoveOutliers_AngleError(sfm_data, 2.0);

  std::cout
    << "\nStructure estimation took (s): " << timer.elapsed() << "." << std::endl
    << "#landmark found: " << sfm_data.GetLandmarks().size() << std::endl;

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
    return EXIT_SUCCESS;
  else
    return EXIT_FAILURE;
}
