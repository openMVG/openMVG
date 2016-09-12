
// Copyright (c) 2015 Pierre MOULON. + modifs by Sébastien CHAPPUIS

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <openMVG/sfm/sfm.hpp>
#include <openMVG/features/features.hpp>
#include <nonFree/sift/SIFT_describer.hpp>
#include <openMVG/image/image.hpp>

#include <openMVG/system/timer.hpp>
#include "openMVG/stl/stl.hpp"

using namespace openMVG;
using namespace openMVG::sfm;

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <cstdlib>

//
// ----------------------------------------------------
//	Image localization from an existing reconstruction
// ----------------------------------------------------
int main(int argc, char **argv)
{
  using namespace std;
  std::cout << std::endl
    << "-----------------------------------------------------------\n"
    << "  Multiple image localization in an existing SfM reconstruction:\n"
    << "-----------------------------------------------------------\n"
    << std::endl;

  CmdLine cmd;

  std::string sSfM_Data_Filename;
  std::string sMatchesDir;
  std::string sOutDir = "";
  std::string sQueryDir;
	bool bSingleIntrinsics = false;
  double dMaxResidualError = std::numeric_limits<double>::infinity();

  cmd.add( make_option('i', sSfM_Data_Filename, "input_file") );
  cmd.add( make_option('m', sMatchesDir, "match_dir") );
  cmd.add( make_option('o', sOutDir, "out_dir") );
  cmd.add( make_option('q', sQueryDir, "query_image_dir"));
  cmd.add( make_option('s', bSingleIntrinsics, "single_intrinsics"));
  cmd.add( make_option('r', dMaxResidualError, "residual_error"));

  try {
    if (argc == 1) throw std::string("Invalid parameter.");
    cmd.process(argc, argv);
  } catch(const std::string& s) {
    std::cerr << "Usage: " << argv[0] << '\n'
    << "[-i|--input_file] path to a SfM_Data scene\n"
    << "[-m|--match_dir] path to the directory containing the matches corresponding to the provided SfM_Data scene\n"
    << "[-o|--out_dir] path where the output data will be stored\n"
    << "[-q|--query_image_dir] path to the directory containing the images that must be localized (can also contain the images from the initial reconstruction)\n"
    << "(optional)\n"
    << "[-r|--residual_error] upper bound of the residual error tolerance\n"
    << "[-s|--single_intrinsics] (default: false) when switched on, the program will check if the input sfm_data contains a single intrinsics and, if so, take\
		this value as input for the intrinsics of the query images.\n"
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

  if ( !stlplus::folder_exists( sQueryDir ) )
  {
    std::cerr << "\nThe query directory doesn't exist" << std::endl;
    return EXIT_FAILURE;
  }

  if (sOutDir.empty())  {
    std::cerr << "\nIt is an invalid output directory" << std::endl;
    return EXIT_FAILURE;
  }

  if (!stlplus::folder_exists(sOutDir))
    stlplus::folder_create(sOutDir);

 if (bSingleIntrinsics && sfm_data.GetIntrinsics().size() != 1)
	{
		std::cout << "More than one intrinsics to compare to in input scene => Consider intrinsics as unkown." << std::endl;
	}

  //-- Localization
  // - init the retrieval database
  // - Go along the sfm_data view
  // - extract the regions of the view
  // - try to locate the images
  // -

  std::vector<Vec3> vec_found_poses;

  sfm::SfM_Localization_Single_3DTrackObservation_Database localizer;
  if (!localizer.Init(sfm_data, *regions_provider.get()))
  {
    std::cerr << "Cannot initialize the SfM localizer" << std::endl;
  }
  // Since we have copied interesting data, release some memory
  regions_provider.reset();

	// list images from sfm_data in vector
	std::vector<std::string> vec_image_original (sfm_data.GetViews().size());
	int n(-1);
	std::generate(vec_image_original.begin(),vec_image_original.end(),[&n,&sfm_data]{ n++; return stlplus::filename_part(sfm_data.views[n].get()->s_Img_path);} );
	
	// list images in query directory
  std::vector<std::string> vec_image = stlplus::folder_files( sQueryDir );
  std::sort(vec_image.begin(), vec_image.end());

	// find difference between two list of images
	std::vector<std::string> vec_image_new;
	std::set_difference(vec_image.begin(),vec_image.end(),vec_image_original.begin(),vec_image_original.end(),std::back_inserter(vec_image_new));

	// copy sfm_data
	SfM_Data sfm_data_out = sfm_data;
	
	// references
	Views & views = sfm_data_out.views;
	Poses & poses = sfm_data_out.poses;
	Intrinsics & intrinsics = sfm_data_out.intrinsics;
	Landmarks & structure = sfm_data_out.structure;
	
	// first intrinsics of the input sfm_data file, to be used if we inforce single intrinsics
	cameras::Pinhole_Intrinsic_Radial_K3 * ptrPinhole = dynamic_cast<cameras::Pinhole_Intrinsic_Radial_K3*>(sfm_data.GetIntrinsics().at(0).get());

  for ( std::vector<std::string>::const_iterator iter_image = vec_image_new.begin();
    iter_image != vec_image_new.end();
    ++iter_image)
  {

    std::unique_ptr<Regions> query_regions(regions_type->EmptyClone());
    image::Image<unsigned char> imageGray;
    {
			string sView_filename = stlplus::create_filespec(sQueryDir,*iter_image);
      if (!image::ReadImage(sView_filename.c_str(), &imageGray))
      {
        std::cerr << "Cannot open the input provided image : " << *iter_image << std::endl;
        return EXIT_FAILURE;
      }
			const std::string
				sFeat = stlplus::create_filespec(sMatchesDir, stlplus::basename_part(sView_filename.c_str()), "feat"),
				sDesc = stlplus::create_filespec(sMatchesDir, stlplus::basename_part(sView_filename.c_str()), "desc");

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

		std::shared_ptr<cameras::IntrinsicBase> optional_intrinsic(nullptr);
		if (bSingleIntrinsics && sfm_data.GetIntrinsics().size() == 1)
		{
			optional_intrinsic = std::make_shared<cameras::Pinhole_Intrinsic_Radial_K3>(
				imageGray.Width(), imageGray.Height(),
				ptrPinhole->focal(), ptrPinhole->principal_point()[0], ptrPinhole->principal_point()[1],
				ptrPinhole->getParams()[3], ptrPinhole->getParams()[4], ptrPinhole->getParams()[5]);
		}

    geometry::Pose3 pose;
    sfm::Image_Localizer_Match_Data matching_data;
    matching_data.error_max = dMaxResidualError;

		// Build the view corresponding to the image
    View v(*iter_image, views.size(), views.size(), views.size(), imageGray.Width(), imageGray.Height());

    // Try to localize the image in the database thanks to its regions
    if (!localizer.Localize(
      Pair(imageGray.Width(), imageGray.Height()),
      optional_intrinsic.get(),
      *(query_regions.get()),
      pose,
      &matching_data))
    {
      std::cerr << "Cannot locate the image " << *iter_image << std::endl;
      v.id_intrinsic = UndefinedIndexT;
      v.id_pose = UndefinedIndexT;
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

			// Add the computed intrinsic to the sfm_container
			intrinsics[v.id_intrinsic] = optional_intrinsic;
			
			// Add the computed pose to the sfm_container
			poses[v.id_pose] = pose;
    }
		
    // Add the view to the sfm_container
    views[v.id_view] = std::make_shared<View>(v);
	}

  GroupSharedIntrinsics(sfm_data_out);
	
	std::cout << " Total poses successfuly computed : " << vec_found_poses.size() << "/" << vec_image_new.size() << endl;

  // Export the found camera position in a ply.
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

	// Export found camera poses along with original camera poses in a new sfm_data file
  if (!Save(
    sfm_data_out,
    stlplus::create_filespec( sOutDir, "sfm_data_expanded.json" ).c_str(),
    ESfM_Data(ALL)))
  {
    return EXIT_FAILURE;
  }
	// export also as ply
  if (!Save(
    sfm_data_out,
    stlplus::create_filespec( sOutDir, "sfm_data_expanded.ply" ).c_str(),
    ESfM_Data(ALL)))
  {
    return EXIT_FAILURE;
  }

  return EXIT_FAILURE;
}
