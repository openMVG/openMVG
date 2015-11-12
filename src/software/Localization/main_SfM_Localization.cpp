
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <openMVG/sfm/sfm.hpp>
#include <openMVG/features/features.hpp>
#include <nonFree/sift/SIFT_describer.hpp>
#include <openMVG/image/image.hpp>

#include <openMVG/system/timer.hpp>

using namespace openMVG;
using namespace openMVG::sfm;

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <cstdlib>

// ---------------------------------------------------------------------------
// Image localization API sample:
// ---------------------------------------------------------------------------
// - Allow to locate an image to an existing SfM_reconstruction
//   if 3D-2D matches are found
// - A demonstration mode (default):
//   - try to locate all the view of the SfM_Data reconstruction
// ---------------------------------------------------------------------------
//
int main(int argc, char **argv)
{
  using namespace std;
  std::cout << std::endl
    << "-----------------------------------------------------------\n"
    << "  Image localization in an existing SfM reconstruction:\n"
    << "-----------------------------------------------------------\n"
    << std::endl;

  CmdLine cmd;

  std::string sSfM_Data_Filename;
  std::string sMatchesDir;
  std::string sOutDir = "";
  std::string sQueryImage;
  double dMaxResidualError = std::numeric_limits<double>::infinity();

  cmd.add( make_option('i', sSfM_Data_Filename, "input_file") );
  cmd.add( make_option('m', sMatchesDir, "match_dir") );
  cmd.add( make_option('o', sOutDir, "out_dir") );
  cmd.add( make_option('q', sQueryImage, "query_image"));
  cmd.add( make_option('r', dMaxResidualError, "residual_error"));

  try {
    if (argc == 1) throw std::string("Invalid parameter.");
    cmd.process(argc, argv);
  } catch(const std::string& s) {
    std::cerr << "Usage: " << argv[0] << '\n'
    << "[-i|--input_file] path to a SfM_Data scene\n"
    << "[-m|--match_dir] path to the matches that corresponds to the provided SfM_Data scene\n"
    << "[-o|--out_dir] path where the output data will be stored\n"
    << "(optional)\n"
    << "[-q|--query_image] path to the image that must be localized\n"
    << "[-r|--residual_error] upper bound of the residual error tolerance\n"
    << std::endl;

    std::cerr << s << std::endl;
    return EXIT_FAILURE;
  }

  // Load input SfM_Data scene
  SfM_Data sfm_data;
  if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(ALL))) {
    std::cerr << std::endl
      << "The input SfM_Data file \""<< sSfM_Data_Filename << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }

  if (sfm_data.GetPoses().empty() || sfm_data.GetLandmarks().empty())
  {
    std::cerr << std::endl
      << "The input SfM_Data file have not 3D content to match with." << std::endl;
    return EXIT_FAILURE;
  }

  // ---------------
  // Initialization
  // ---------------

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

  // Init the feature extractor that have been used for the reconstruction
  using namespace openMVG::features;
  std::unique_ptr<Image_describer> image_describer;
  if (stlplus::is_file(sImage_describer))
  {
    // Dynamically load the image_describer from the file (will restore old used settings)
    std::ifstream stream(sImage_describer.c_str());
    if (!stream.is_open())
      return false;

    try
    {
      cereal::JSONInputArchive archive(stream);
      archive(cereal::make_nvp("image_describer", image_describer));
    }
    catch (const cereal::Exception & e)
    {
      std::cerr << e.what() << std::endl
        << "Cannot dynamically allocate the Image_describer interface." << std::endl;
      return EXIT_FAILURE;
    }
  }
  else
  {
    std::cerr << "Expected file image_describer.json cannot be opened." << std::endl;
    return EXIT_FAILURE;
  }

  // Load the SfM_Data region's views
  std::shared_ptr<Regions_Provider> regions_provider = std::make_shared<Regions_Provider>();
  if (!regions_provider->load(sfm_data, sMatchesDir, regions_type)) {
    std::cerr << std::endl << "Invalid regions." << std::endl;
    return EXIT_FAILURE;
  }

  if (sOutDir.empty())  {
    std::cerr << "\nIt is an invalid output directory" << std::endl;
    return EXIT_FAILURE;
  }

  if (!stlplus::folder_exists(sOutDir))
    stlplus::folder_create(sOutDir);

  //-
  //-- Localization
  // - init the retrieval database
  // - Go along the sfm_data view
  // - extract the regions of the view
  // - try to locate the image
  //-

  std::vector<Vec3> vec_found_poses;

  sfm::SfM_Localization_Single_3DTrackObservation_Database localizer;
  if (!localizer.Init(sfm_data, *regions_provider.get()))
  {
    std::cerr << "Cannot initialize the SfM localizer" << std::endl;
  }
  // Since we have copied interesting data, release some memory
  regions_provider.reset();

  if (!sQueryImage.empty())
  {
    std::cout << "SfM::localization => try with image: " << sQueryImage << std::endl;
    std::unique_ptr<Regions> query_regions;
    image::Image<unsigned char> imageGray;
    {
      if (!image::ReadImage(sQueryImage.c_str(), &imageGray))
      {
        std::cerr << "Cannot open the input provided image" << std::endl;
        return EXIT_FAILURE;
      }
      // Compute features and descriptors
      image_describer->Describe(imageGray, query_regions);
      std::cout << "#regions detected in query image: " << query_regions->RegionCount() << std::endl;
    }

    // Suppose intrinsic as unkown
    std::shared_ptr<cameras::IntrinsicBase> optional_intrinsic (nullptr);

    geometry::Pose3 pose;
    sfm::Image_Localizer_Match_Data matching_data;
    matching_data.error_max = dMaxResidualError;

    // Try to localize the image in the database thanks to its regions
    if (!localizer.Localize(
      Pair(imageGray.Width(), imageGray.Height()),
      optional_intrinsic.get(),
      *(query_regions.get()),
      pose,
      &matching_data))
    {
      std::cerr << "Cannot locate the image" << std::endl;
    }
    else
    {
      const bool b_new_intrinsic = (optional_intrinsic == nullptr);
      // A valid pose has been found (try to refine it):
      // If not intrinsic as input:
      //  init a new one from the projection matrix decomposition
      // Else use the existing one and consider as static.
      if (b_new_intrinsic)
      {
        // setup a default camera model from the found projection matrix
        Mat3 K, R;
        Vec3 t;
        KRt_From_P(matching_data.projection_matrix, &K, &R, &t);

        const double focal = (K(0,0) + K(1,1))/2.0;
        const Vec2 principal_point(K(0,2), K(1,2));
        optional_intrinsic = std::make_shared<cameras::Pinhole_Intrinsic_Radial_K3>(
          imageGray.Width(), imageGray.Height(),
          focal, principal_point(0), principal_point(1));
      }
      sfm::SfM_Localizer::RefinePose
      (
        optional_intrinsic.get(),
        pose, matching_data,
        true, b_new_intrinsic
      );

      vec_found_poses.push_back(pose.center());
    }
  }
  else
  {
    std::cout << "----------------\n"
      << " DEMONSTRATION \n"
      << "----------------\n" << std::endl;
    std::cout << "Will try to locate all the view of the loaded sfm_data file." << std::endl;

    for (const auto & viewIter : sfm_data.GetViews())
    {
      const View * view = viewIter.second.get();
      // Load an image, extract the regions and match
      const std::string sImagePath = sfm_data.s_root_path + view->s_Img_path;

      std::cout << "SfM::localization => try with image: " << sImagePath << std::endl;

      std::unique_ptr<Regions> query_regions(regions_type->EmptyClone());
      const std::string basename = stlplus::basename_part(sImagePath);
      const std::string featFile = stlplus::create_filespec(sMatchesDir, basename, ".feat");
      const std::string descFile = stlplus::create_filespec(sMatchesDir, basename, ".desc");
      if (!query_regions->Load(featFile, descFile))
      {
        std::cerr << "Invalid regions files for the view: " << sImagePath << std::endl;
        continue;
      }

      // Initialize intrinsics data for the view if any
      std::shared_ptr<cameras::IntrinsicBase> optional_intrinsic (nullptr);
      if (sfm_data.GetIntrinsics().count(view->id_intrinsic))
      {
        optional_intrinsic = sfm_data.GetIntrinsics().at(view->id_intrinsic);
      }

      geometry::Pose3 pose;
      sfm::Image_Localizer_Match_Data matching_data;
      matching_data.error_max = dMaxResidualError;

      // Try to localize the image in the database thanks to its regions
      if (!localizer.Localize(
        Pair(view->ui_width, view->ui_height),
        optional_intrinsic.get(),
        *(query_regions.get()),
        pose,
        &matching_data))
      {
        std::cerr << "Cannot locate the image" << std::endl;
      }
      else
      {
        const bool b_new_intrinsic = (optional_intrinsic == nullptr);
        // A valid pose has been found (try to refine it):
        // If no valid intrinsic as input:
        //  init a new one from the projection matrix decomposition
        // Else use the existing one and consider it as constant
        if (b_new_intrinsic)
        {
          // setup a default camera model from the found projection matrix
          Mat3 K, R;
          Vec3 t;
          KRt_From_P(matching_data.projection_matrix, &K, &R, &t);

          const double focal = (K(0,0) + K(1,1))/2.0;
          const Vec2 principal_point(K(0,2), K(1,2));
          optional_intrinsic = std::make_shared<cameras::Pinhole_Intrinsic_Radial_K3>(
            view->ui_width, view->ui_height,
            focal, principal_point(0), principal_point(1));
        }
        sfm::SfM_Localizer::RefinePose
        (
          optional_intrinsic.get(),
          pose, matching_data,
          true, b_new_intrinsic
        );

        vec_found_poses.push_back(pose.center());
      }
    }
    std::cout
      << "\n#images found: " << vec_found_poses.size()
      << "\n#sfm_data view count: " << sfm_data.GetViews().size()
      << std::endl;
  }

  // Export the found camera position
  const std::string out_file_name = stlplus::create_filespec(sOutDir, "found_pose_centers", "ply");
  {
    std::ofstream outfile;
    outfile.open(out_file_name.c_str(), std::ios_base::out);
    if (outfile.is_open()) {
      outfile << "ply"
       << "\n" << "format ascii 1.0"
       << "\n" << "element vertex " << vec_found_poses.size()
       << "\n" << "property float x"
       << "\n" << "property float y"
       << "\n" << "property float z"
       << "\n" << "property uchar red"
       << "\n" << "property uchar green"
       << "\n" << "property uchar blue"
       << "\n" << "end_header" << "\n";

      for (const Vec3 & pose_center: vec_found_poses) {
        outfile << pose_center.transpose() << " " << "255 0 0" << "\n";
      }

      outfile.close();
    }
  }

  return EXIT_FAILURE;
}
