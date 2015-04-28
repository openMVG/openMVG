
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include "openMVG/sfm/sfm.hpp"
#include "openMVG/system/timer.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

using namespace openMVG;

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

  cmd.add( make_option('i', sSfM_Data_Filename, "input_file") );
  cmd.add( make_option('m', sMatchesDir, "match_dir") );
  cmd.add( make_option('f', sMatchFile, "match_file") );
  cmd.add( make_option('o', sOutFile, "output_file") );

  try {
    if (argc == 1) throw std::string("Invalid command line parameter.");
    cmd.process(argc, argv);
  } catch(const std::string& s) {
    std::cerr << "Usage: " << argv[0] << '\n'
    << "[-i|--input_file] path to a SfM_Data scene\n"
    << "[-m|--match_dir] path to the features and descriptor that "
    << " corresponds to the provided SfM_Data scene\n"
    << "[-f|--match_file] (opt.) path to a matches file (used pairs will be used)\n"
    << "[-o|--output_file] path where the output data will be stored\n"
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

  // Prepare the features provider
  std::shared_ptr<Features_Provider> feats_provider = std::make_shared<Features_Provider>();
  if (!feats_provider->load(sfm_data, sMatchesDir)) {
    std::cerr << std::endl
      << "Invalid features." << std::endl;
    return EXIT_FAILURE;
  }

  //--
  //- Pair selection method
  //  - geometry guided -> camera frustum intersection
  //  - putative matches guided (photometric matches) (Keep pair that have valid Intrinsic & Pose ids)
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
    if (!matching::PairedIndMatchImport(sMatchFile, matches)) {
      std::cerr<< "Unable to read the matches file." << std::endl;
      return EXIT_FAILURE;
    }
    pairs = getPairs(matches);
  }

  // Keep only Pairs that belong to valid view indexes.
  std::set<IndexT> valid_viewIdx = Get_Valid_Views(sfm_data);
  pairs = Pair_filter(pairs, valid_viewIdx);

  // Collect view features descriptors.
  typedef Descriptor<unsigned char, 128> DescriptorT;
  typedef vector<DescriptorT > DescsT;

  Hash_Map<IndexT, DescsT > desc_per_view;
  for(Views::const_iterator iterViews = sfm_data.views.begin();
      iterViews != sfm_data.views.end();
      ++iterViews)
  {
    const View * view = iterViews->second.get();
    const std::string sView_filename = stlplus::create_filespec(sfm_data.s_root_path,
      view->s_Img_path);
    const std::string sDescJ = stlplus::create_filespec(sMatchesDir,
      stlplus::basename_part(sView_filename), "desc");

    loadDescsFromBinFile(sDescJ, desc_per_view[view->id_view]);
  }

  openMVG::Timer timer;

  //------------------------------------------
  // Compute Structure from known camera poses
  //------------------------------------------
  SfM_Data_Structure_Estimation_From_Known_Poses structure_estimator;
  structure_estimator.run(sfm_data, pairs, feats_provider.get(), desc_per_view);
  RemoveOutliers_AngleError(sfm_data, 2.0);

  std::cout << "\nStructure estimation took (s): " << timer.elapsed() << "." << std::endl;

  std::cout << "#landmark found: " << sfm_data.getLandmarks().size() << std::endl;

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
