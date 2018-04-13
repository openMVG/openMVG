// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// The <cereal/archives> headers are special and must be included first.
#include <cereal/archives/json.hpp>

#include <openMVG/sfm/sfm.hpp>
#include <openMVG/features/feature.hpp>
#include <openMVG/features/image_describer.hpp>
#include <openMVG/image/image_io.hpp>
#include <software/SfM/SfMPlyHelper.hpp>

#include <openMVG/system/timer.hpp>
#include "openMVG/stl/stl.hpp"

using namespace openMVG;
using namespace openMVG::sfm;

#include "nonFree/sift/SIFT_describer_io.hpp"
#include "openMVG/features/akaze/image_describer_akaze_io.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <cstdlib>

#ifdef OPENMVG_USE_OPENMP
#include <omp.h>
#endif

// Naive function for finding the biggest common root dir from two paths
std::string FindCommonRootDir(const std::string & dir1, const std::string & dir2)
{
  int i = 0;
  for (; i != std::min(dir1.size(), dir2.size()); i++)
  {
    if (dir1[i] != dir2[i]) break;
  }
  return dir1.substr(0,i);
}

// ----------------------------------------------------
// Multiple Images localization from an existing reconstruction
// ----------------------------------------------------
int main(int argc, char **argv)
{
  using namespace std;
  std::cout << std::endl
    << "-----------------------------------------------------------\n"
    << "  Images localization in an existing SfM reconstruction:\n"
    << "-----------------------------------------------------------\n"
    << std::endl;

  CmdLine cmd;

  std::string sSfM_Data_Filename;
  std::string sMatchesDir;
  std::string sOutDir = "";
  std::string sMatchesOutDir;
  std::string sQueryDir;
  double dMaxResidualError = std::numeric_limits<double>::infinity();
  int i_User_camera_model = cameras::PINHOLE_CAMERA_RADIAL3;
  bool bUseSingleIntrinsics = false;
  bool bExportStructure = false;

#ifdef OPENMVG_USE_OPENMP
  int iNumThreads = 0;
#endif

  cmd.add( make_option('i', sSfM_Data_Filename, "input_file") );
  cmd.add( make_option('m', sMatchesDir, "match_dir") );
  cmd.add( make_option('o', sOutDir, "out_dir") );
  cmd.add( make_option('u', sMatchesOutDir, "match_out_dir") );
  cmd.add( make_option('q', sQueryDir, "query_image_dir"));
  cmd.add( make_option('r', dMaxResidualError, "residual_error"));
  cmd.add( make_option('c', i_User_camera_model, "camera_model") );
  cmd.add( make_switch('s', "single_intrinsics"));
  cmd.add( make_switch('e', "export_structure"));

#ifdef OPENMVG_USE_OPENMP
  cmd.add( make_option('n', iNumThreads, "numThreads") );
#endif


  try {
    if (argc == 1) throw std::string("Invalid parameter.");
    cmd.process(argc, argv);
  } catch (const std::string& s) {
    std::cerr << "Usage: " << argv[0] << '\n'
    << "[-i|--input_file] path to a SfM_Data scene\n"
    << "[-m|--match_dir] path to the directory containing the matches\n"
    << "  corresponding to the provided SfM_Data scene\n"
    << "[-o|--out_dir] path where the output data will be stored\n"
    << "[-u|--match_out_dir] path to the directory where new matches will be stored\n"
    << "  (if empty the initial matching directory will be used)\n"
    << "[-q|--query_image_dir] path to an image OR to the directory containing the images that must be localized\n"
    << "  (the directory can also contain the images from the initial reconstruction)\n"
    << "\n"
    << "(optional)\n"
    << "[-r|--residual_error] upper bound of the residual error tolerance\n"
    << "[-s|--single_intrinsics] (switch) when switched on, the program will check if the input sfm_data\n"
    << "  contains a single intrinsics and, if so, take this value as intrinsics for the query images.\n"
    << "  (OFF by default)\n"
    << "[-e|--export_structure] (switch) when switched on, the program will also export structure to output sfm_data.\n"
    << "  if OFF only VIEWS, INTRINSICS and EXTRINSICS are exported (OFF by default)\n"
    << "[-c|--camera_model] Camera model type for view with unknown intrinsic:\n"
      << "\t 1: Pinhole\n"
      << "\t 2: Pinhole radial 1\n"
      << "\t 3: Pinhole radial 3 (default)\n"
      << "\t 4: Pinhole radial 3 + tangential 2\n"
      << "\t 5: Pinhole fisheye\n"
      << "\t 7: Spherical camera\n"
#ifdef OPENMVG_USE_OPENMP
    << "[-n|--numThreads] number of thread(s)\n"
#endif
    << std::endl;

    std::cerr << s << std::endl;
    return EXIT_FAILURE;
  }

  if ( !isValid(openMVG::cameras::EINTRINSIC(i_User_camera_model)) )  {
    std::cerr << "\n Invalid camera type" << std::endl;
    return EXIT_FAILURE;
  }

  // Load input SfM_Data scene
  SfM_Data sfm_data;
  if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(ALL))) {
    std::cerr << std::endl
      << "The input SfM_Data file \""<< sSfM_Data_Filename << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }

  if (sMatchesOutDir.empty())
  {
    sMatchesOutDir = sMatchesDir;
  }

  if (sfm_data.GetPoses().empty() || sfm_data.GetLandmarks().empty())
  {
    std::cerr << std::endl
      << "The input SfM_Data file have not 3D content to match with." << std::endl;
    return EXIT_FAILURE;
  }

  bUseSingleIntrinsics = cmd.used('s');
  bExportStructure = cmd.used('e');
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
  std::unique_ptr<Image_describer> image_describer;
  if (stlplus::is_file(sImage_describer))
  {
    // Dynamically load the image_describer from the file (will restore old used settings)
    std::ifstream stream(sImage_describer.c_str());
    if (!stream.is_open())
      return EXIT_FAILURE;

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

  // Show the progress on the command line:
  C_Progress_display progress;

  // Load the SfM_Data region's views
  std::shared_ptr<Regions_Provider> regions_provider = std::make_shared<Regions_Provider>();
  if (!regions_provider->load(sfm_data, sMatchesDir, regions_type, &progress)) {
    std::cerr << std::endl << "Invalid regions." << std::endl;
    return EXIT_FAILURE;
  }

  if ( !stlplus::folder_exists( sQueryDir ) && !stlplus::file_exists( sQueryDir ) )
  {
    std::cerr << "\nThe query directory/file does not exist : " << std::endl;
    std::cerr << sQueryDir << std::endl;
    return EXIT_FAILURE;
  }

  if (sOutDir.empty())  {
    std::cerr << "\nPlease provide a valid directory for the option [-o|--out_dir]." << std::endl;
    return EXIT_FAILURE;
  }

  if (!stlplus::folder_exists(sOutDir))
    stlplus::folder_create(sOutDir);

  if (bUseSingleIntrinsics && sfm_data.GetIntrinsics().size() != 1)
  {
    std::cout << "More than one intrinsics to compare to in input scene "
              << " => Consider intrinsics as unkown." << std::endl;
  }

  //-- Localization
  // - init the retrieval database
  // - Go along the sfm_data view
  // - extract the regions of the view
  // - try to locate the images
  // - add the images to the sfm_data scene

  std::vector<Vec3> vec_found_poses;

  sfm::SfM_Localization_Single_3DTrackObservation_Database localizer;
  if (!localizer.Init(sfm_data, *regions_provider.get()))
  {
    std::cerr << "Cannot initialize the SfM localizer" << std::endl;
  }
  // Since we have copied interesting data, release some memory
  regions_provider.reset();

  // list images from sfm_data in a vector
  std::vector<std::string> vec_image_original (sfm_data.GetViews().size());
  int n(-1);
  std::generate(vec_image_original.begin(),
                vec_image_original.end(),
                [&n,&sfm_data]
                {
                  n++;
                  return stlplus::filename_part(sfm_data.views.at(n)->s_Img_path);
                });

  // list images in query directory
  std::vector<std::string> vec_image;

  if (stlplus::is_file(sQueryDir))
  {
    vec_image.emplace_back(stlplus::filename_part(sQueryDir)); // single file
    sQueryDir = stlplus::folder_part(sQueryDir);
  }
  else vec_image = stlplus::folder_files(sQueryDir); // multiple files

  std::sort(vec_image.begin(), vec_image.end());

  // find difference between two list of images
  std::vector<std::string> vec_image_new;
  std::set_difference(vec_image.cbegin(), vec_image.cend(),
      vec_image_original.cbegin(),vec_image_original.cend(),
      std::back_inserter(vec_image_new));

  // find common root directory between images in vec_image_original and vec_images_new
  const std::string common_root_dir = FindCommonRootDir(sfm_data.s_root_path, sQueryDir);

  // check if sfm_data's root dir differs from the common root dir.
  if (sfm_data.s_root_path != common_root_dir)
  {
    // in that case we have to change all the image paths from the original
    // reconstruction
    for (auto & view : sfm_data.GetViews())
    {
      view.second->s_Img_path = stlplus::create_filespec(stlplus::folder_to_relative_path(common_root_dir, sfm_data.s_root_path),
          view.second->s_Img_path);
    }
    // change root path to common root path
    sfm_data.s_root_path = common_root_dir;
  }

  // references
  Views & views = sfm_data.views;
  Poses & poses = sfm_data.poses;
  Intrinsics & intrinsics = sfm_data.intrinsics;

  int total_num_images = 0;

#ifdef OPENMVG_USE_OPENMP
  const unsigned int nb_max_thread = (iNumThreads == 0) ? 0 : omp_get_max_threads();
    omp_set_num_threads(nb_max_thread);
    #pragma omp parallel for schedule(dynamic)
#endif
  for (int i = 0; i < static_cast<int>(vec_image_new.size()); ++i)
  {
    std::vector<std::string>::const_iterator iter_image = vec_image_new.begin();
    std::advance(iter_image, i);


    // Test if the image format is supported:
    if (openMVG::image::GetFormat((*iter_image).c_str()) == openMVG::image::Unknown)
    {
      std::cerr << *iter_image << " : unknown image file format." << std::endl;
      continue;
    }

    std::cout << "SfM::localization => try with image: " << *iter_image << std::endl;
    std::unique_ptr<Regions> query_regions(regions_type->EmptyClone());
    image::Image<unsigned char> imageGray;
    {
      const std::string sView_filename = stlplus::create_filespec(sQueryDir, *iter_image);
      // Try to open image
      if (!image::ReadImage(sView_filename.c_str(), &imageGray))
      {
        std::cerr << "Cannot open the input provided image : " << *iter_image << std::endl;
        continue;
      }

      const std::string
        sFeat = stlplus::create_filespec(sMatchesOutDir, stlplus::basename_part(sView_filename.c_str()), "feat"),
        sDesc = stlplus::create_filespec(sMatchesOutDir, stlplus::basename_part(sView_filename.c_str()), "desc");

      // Compute features and descriptors and save them if they don't exist yet
      if (!stlplus::file_exists(sFeat) || !stlplus::file_exists(sDesc))
      {
        image_describer->Describe(imageGray, query_regions);
        image_describer->Save(query_regions.get(), sFeat, sDesc);
        std::cout << "#regions detected in query image: " << query_regions->RegionCount() << std::endl;
      }
      else // load already existing regions
      {
        query_regions->Load(sFeat,sDesc);
      }
    }

    std::shared_ptr<cameras::IntrinsicBase> optional_intrinsic;
    if (bUseSingleIntrinsics)
    {
      if (sfm_data.GetIntrinsics().size() != 1)
      {
        std::cerr << "You choose the single intrinsic mode but the sfm_data scene,"
          <<" have too few or too much intrinsics."
          << std::endl;
        continue;
      }
      optional_intrinsic = sfm_data.GetIntrinsics().at(0);
      if (imageGray.Width() != optional_intrinsic->w() || optional_intrinsic->h() != imageGray.Height())
      {
        std::cout << "The provided image does not have the same size as the camera model you want to use." << std::endl;
        continue;
      }
    }
    if (optional_intrinsic)
    {
      std::cout << "- use known intrinsics." << std::endl;
    }
    else
    {
      std::cout << "- use Unknown intrinsics for the resection. A new camera (intrinsic) will be created." << std::endl;

      // Since the spherical image is only defined by its image size we can initialize its camera model.
      // This way the resection will be performed with valid bearing vector
      if (openMVG::cameras::EINTRINSIC(i_User_camera_model) == cameras::CAMERA_SPHERICAL)
      {
        optional_intrinsic = std::make_shared<cameras::Intrinsic_Spherical>(imageGray.Width(), imageGray.Height());
      }
    }

    geometry::Pose3 pose;
    sfm::Image_Localizer_Match_Data matching_data;
    matching_data.error_max = dMaxResidualError;

    bool bSuccessfulLocalization = false;

    // Try to localize the image in the database thanks to its regions
    if (!localizer.Localize(
      optional_intrinsic ? resection::SolverType::P3P_KE_CVPR17 : resection::SolverType::DLT_6POINTS,
      {imageGray.Width(), imageGray.Height()},
      optional_intrinsic.get(),
      *(query_regions.get()),
      pose,
      &matching_data))
    {
      std::cerr << "Cannot locate the image " << *iter_image << std::endl;
      bSuccessfulLocalization = false;
    }
    else
    {
      const bool b_new_intrinsic = (optional_intrinsic == nullptr);
      // A valid pose has been found (try to refine it):
      // If not intrinsic as input:
      // init a new one from the projection matrix decomposition
      // Else use the existing one and consider as static.
      if (b_new_intrinsic)
      {
        // setup a default camera model from the found projection matrix
        Mat3 K, R;
        Vec3 t;
        KRt_From_P(matching_data.projection_matrix, &K, &R, &t);

        const double focal = (K(0,0) + K(1,1))/2.0;
        const Vec2 principal_point(K(0,2), K(1,2));

        switch (openMVG::cameras::EINTRINSIC(i_User_camera_model))
        {
          case cameras::PINHOLE_CAMERA:
            optional_intrinsic = std::make_shared<cameras::Pinhole_Intrinsic>(imageGray.Width(), imageGray.Height(),focal, principal_point(0), principal_point(1));
          break;
          case cameras::PINHOLE_CAMERA_RADIAL1:
            optional_intrinsic = std::make_shared<cameras::Pinhole_Intrinsic_Radial_K1>(imageGray.Width(), imageGray.Height(),focal, principal_point(0), principal_point(1));
          break;
          case cameras::PINHOLE_CAMERA_RADIAL3:
            optional_intrinsic = std::make_shared<cameras::Pinhole_Intrinsic_Radial_K3>(imageGray.Width(), imageGray.Height(),focal, principal_point(0), principal_point(1));
          break;
          case cameras::PINHOLE_CAMERA_BROWN:
            optional_intrinsic = std::make_shared<cameras::Pinhole_Intrinsic_Brown_T2>(imageGray.Width(), imageGray.Height(),focal, principal_point(0), principal_point(1));
          break;
          case cameras::PINHOLE_CAMERA_FISHEYE:
            optional_intrinsic = std::make_shared<cameras::Pinhole_Intrinsic_Fisheye>(imageGray.Width(), imageGray.Height(),focal, principal_point(0), principal_point(1));
          break;
          case cameras::CAMERA_SPHERICAL:
            std::cerr << "The spherical camera cannot be created there. Resection of a spherical camera must be done with an existing camera model." << std::endl;
          break;
          default:
            std::cerr << "Error: unknown camera model: " << static_cast<int>(i_User_camera_model) << std::endl;
        }
      }
      if (optional_intrinsic && sfm::SfM_Localizer::RefinePose(
        optional_intrinsic.get(),
        pose, matching_data,
        true, b_new_intrinsic))
      {
        bSuccessfulLocalization = true;
      }
      else
      {
        std::cerr << "Refining pose for the image " << *iter_image << " failed." << std::endl;
      }

    }
#ifdef OPENMVG_USE_OPENMP
    #pragma omp critical
#endif
    {
      total_num_images++;

      View v(*iter_image, views.size(), views.size(), views.size(), imageGray.Width(), imageGray.Height());
      if (bSuccessfulLocalization)
      {
        vec_found_poses.push_back(pose.center());

        // Add the computed intrinsic to the sfm_container
        if (!bUseSingleIntrinsics)
          intrinsics[v.id_intrinsic] = optional_intrinsic;
        else // Make the view using the existing intrinsic id
          v.id_intrinsic = sfm_data.GetViews().begin()->second->id_intrinsic;
        // Add the computed pose to the sfm_container
        poses[v.id_pose] = pose;

      }
      else
      {
        v.id_intrinsic = UndefinedIndexT;
        v.id_pose = UndefinedIndexT;
      }
      // Add the view to the sfm_container
      views[v.id_view] = std::make_shared<View>(v);
    }
  }

  GroupSharedIntrinsics(sfm_data);

  std::cout << " Total poses found : " << vec_found_poses.size() << "/" << total_num_images << endl;

  // Export the found camera position in a ply.
  const std::string out_file_name = stlplus::create_filespec(sOutDir, "found_pose_centers", "ply");
  plyHelper::exportToPly(vec_found_poses, out_file_name);

  // Export found camera poses along with original reconstruction in a new sfm_data file
  ESfM_Data flag_save;
  if (bExportStructure)
  {
    flag_save = ESfM_Data(ALL);
  }
  else
  {
    flag_save = ESfM_Data(VIEWS|INTRINSICS|EXTRINSICS);
  }
  if (!Save(
    sfm_data,
    stlplus::create_filespec( sOutDir, "sfm_data_expanded.json" ).c_str(),
    flag_save))
  {
    return EXIT_FAILURE;
  }
  // export also as ply
  if (!Save(
    sfm_data,
    stlplus::create_filespec( sOutDir, "sfm_data_expanded.ply" ).c_str(),
    ESfM_Data(ALL)))
  {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
