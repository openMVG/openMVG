
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/image/image.hpp"
#include "openMVG/sfm/sfm.hpp"

/// Feature/Regions & Image describer interfaces
#include "openMVG/features/features.hpp"
#include "nonFree/sift/SIFT_describer.hpp"
#include "nonFree/sift/SIFT_float_describer.hpp"

#if HAVE_CCTAG
#include "openMVG/features/cctag/CCTAG_describer.hpp"
#include "openMVG/features/cctag/SIFT_CCTAG_describer.hpp"
#endif

#include <cereal/archives/json.hpp>
#include "openMVG/system/timer.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "third_party/progress/progress.hpp"

#include <cstdlib>
#include <fstream>
#include <iostream>



using namespace openMVG;
using namespace openMVG::image;
using namespace openMVG::features;
using namespace openMVG::sfm;
using namespace std;


/// - Compute view image description (feature & descriptor extraction)
/// - Export computed data
int main(int argc, char **argv)
{
  CmdLine cmd;

  std::string sSfM_Data_Filename;
  std::string sOutDir = "";
  bool bUpRight = false;
  std::string sImage_Describer_Method = "SIFT";
  bool bForce = false;
  std::string sFeaturePreset = "";
  int rangeStart = -1;
  int rangeSize = 1;

  // required
  cmd.add( make_option('i', sSfM_Data_Filename, "input_file") );
  cmd.add( make_option('o', sOutDir, "outdir") );
  // Optional
  cmd.add( make_option('m', sImage_Describer_Method, "describerMethod") );
  cmd.add( make_option('u', bUpRight, "upright") );
  cmd.add( make_option('f', bForce, "force") );
  cmd.add( make_option('p', sFeaturePreset, "describerPreset") );
  cmd.add( make_option('s', rangeStart, "range_start") );
  cmd.add( make_option('r', rangeSize, "range_size") );

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch(const std::string& s) {
      std::cerr << "Usage: " << argv[0] << '\n'
      << "[-i|--input_file] a SfM_Data file or a directory with JPG files\n"
      << "[-o|--outdir path] \n"
      << "\n[Optional]\n"
      << "[-f|--force] Force to recompute data\n"
      << "[-m|--describerMethod]\n"
      << "  (method to use to describe an image):\n"
      << "   SIFT (default),\n"
      << "   SIFT_FLOAT to use SIFT stored as float,\n"
      << "   AKAZE_FLOAT: AKAZE with floating point descriptors,\n"
      << "   AKAZE_MLDB:  AKAZE with binary descriptors]\n"
#if HAVE_CCTAG
      << "   CCTAG3: CCTAG markers with 3 crowns\n"
      << "   CCTAG3: CCTAG markers with 4 crowns\n"
      << "   SIFT_CCTAG3: CCTAG markers with 3 crowns\n" 
      << "   SIFT_CCTAG4: CCTAG markers with 4 crowns\n" 
#endif
      << "[-u|--upright] Use Upright feature 0 or 1\n"
      << "[-p|--describerPreset]\n"
      << "  (used to control the Image_describer configuration):\n"
      << "   LOW,\n"
      << "   MEDIUM,\n"
      << "   NORMAL (default),\n"
      << "   HIGH,\n"
      << "   ULTRA: !!Can take long time!!\n"
      << "[-s]--range_start] range image index start\n"
      << "[-r]--range_size] range size\n"
      << std::endl;

      std::cerr << s << std::endl;
      return EXIT_FAILURE;
  }

  std::cout << " You called : " <<std::endl
            << argv[0] << std::endl
            << "--input_file " << sSfM_Data_Filename << std::endl
            << "--outdir " << sOutDir << std::endl
            << "--describerMethod " << sImage_Describer_Method << std::endl
            << "--upright " << bUpRight << std::endl
            << "--describerPreset " << (sFeaturePreset.empty() ? "NORMAL" : sFeaturePreset) << std::endl
            << "--force " << bForce << std::endl
            << "--range_start " << rangeStart << std::endl
            << "--range_size " << rangeSize << std::endl;


  if (sOutDir.empty())  {
    std::cerr << "\nIt is an invalid output directory" << std::endl;
    return EXIT_FAILURE;
  }

  // Create output dir
  if (!stlplus::folder_exists(sOutDir))
  {
    if (!stlplus::folder_create(sOutDir))
    {
      std::cerr << "Cannot create output directory" << std::endl;
      return EXIT_FAILURE;
    }
  }

  //---------------------------------------
  // a. Load input scene
  //---------------------------------------
  SfM_Data sfm_data;

  if(sSfM_Data_Filename.empty())
  {
    std::cerr << "\nThe input directory doesn't exist" << std::endl;
    return EXIT_FAILURE;
  }
  else if(stlplus::is_file( sSfM_Data_Filename))
  {
    if(!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(VIEWS|INTRINSICS)))
    {
      std::cerr << std::endl
        << "The input file \""<< sSfM_Data_Filename << "\" cannot be read" << std::endl;
      return EXIT_FAILURE;
    }
  }
  else if(stlplus::is_folder(sSfM_Data_Filename))
  {
    // Retrieve image paths
    std::vector<std::string> vec_images;
    const std::vector<std::string> supportedExtensions {"jpg", "jpeg"};
    
    vec_images = stlplus::folder_files(sSfM_Data_Filename);
    std::sort(vec_images.begin(), vec_images.end());
    
    sfm_data.s_root_path = "";
    if(!sSfM_Data_Filename.empty())
      sfm_data.s_root_path = sSfM_Data_Filename; // Setup main image root_path
    Views & views = sfm_data.views;
    
    for(std::vector<std::string>::const_iterator iter_image = vec_images.begin();
        iter_image != vec_images.end(); ++iter_image)
    {
      std::size_t id_view = 0;
      stl::hash_combine(id_view, *iter_image);
    
      // Build the view corresponding to the image
      View v(*iter_image, (IndexT)id_view);
      v.id_intrinsic = UndefinedIndexT;
      views[v.id_view] = std::make_shared<View>(v);
    }
  }

  // b. Init the image_describer
  // - retrieve the used one in case of pre-computed features
  // - else create the desired one

  using namespace openMVG::features;
  std::unique_ptr<Image_describer> image_describer;

  const std::string sImage_describer = stlplus::create_filespec(sOutDir, "image_describer", "json");
  if (!bForce && stlplus::is_file(sImage_describer))
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
    // Create the desired Image_describer method.
    // Don't use a factory, perform direct allocation
    if (sImage_Describer_Method == "SIFT")
    {
      image_describer.reset(new SIFT_Image_describer(SiftParams(), !bUpRight));
    }
    else
    if (sImage_Describer_Method == "SIFT_FLOAT")
    {
      image_describer.reset(new SIFT_float_describer(SiftParams(), !bUpRight));
    }
#if HAVE_CCTAG
    else
    if (sImage_Describer_Method == "CCTAG3")
    {
      image_describer.reset(new CCTAG_Image_describer(3));
    }
    else
    if (sImage_Describer_Method == "CCTAG4")
    {
      image_describer.reset(new CCTAG_Image_describer(4));
    }
    else
    if (sImage_Describer_Method == "SIFT_CCTAG3")
    {
      image_describer.reset(new SIFT_CCTAG_Image_describer(SiftParams(), !bUpRight, 3));
    }
    else
    if (sImage_Describer_Method == "SIFT_CCTAG4")
    {
      image_describer.reset(new SIFT_CCTAG_Image_describer(SiftParams(), !bUpRight, 4));
    }
#endif //HAVE_CCTAG   
    else
    if (sImage_Describer_Method == "AKAZE_FLOAT")
    {
      image_describer.reset(new AKAZE_Image_describer(AKAZEParams(AKAZEConfig(), AKAZE_MSURF), !bUpRight));
    }
    else
    if (sImage_Describer_Method == "AKAZE_MLDB")
    {
      image_describer.reset(new AKAZE_Image_describer(AKAZEParams(AKAZEConfig(), AKAZE_MLDB), !bUpRight));
    }
    //image_describer.reset(new AKAZE_Image_describer(AKAZEParams(AKAZEConfig(), AKAZE_LIOP), !bUpRight));
    if (!image_describer)
    {
      std::cerr << "Cannot create the designed Image_describer:"
        << sImage_Describer_Method << "." << std::endl;
      return EXIT_FAILURE;
    }
    else
    {
      if (!sFeaturePreset.empty())
      if (!image_describer->Set_configuration_preset(sFeaturePreset))
      {
        std::cerr << "Preset configuration failed." << std::endl;
        return EXIT_FAILURE;
      }
    }

    // Export the used Image_describer and region type for:
    // - dynamic future regions computation and/or loading
    {
      std::ofstream stream(sImage_describer.c_str());
      if (!stream.is_open())
        return false;

      cereal::JSONOutputArchive archive(stream);
      archive(cereal::make_nvp("image_describer", image_describer));
      std::unique_ptr<Regions> regionsType;
      image_describer->Allocate(regionsType);
      archive(cereal::make_nvp("regions_type", regionsType));
    }
  }

  // Feature extraction routines
  // For each View of the SfM_Data container:
  // - if regions file exist continue,
  // - if no file, compute features
  {
    system::Timer timer;
    C_Progress_display my_progress_bar( sfm_data.GetViews().size(),
      std::cout, "\n- EXTRACT FEATURES -\n" );

    Views::const_iterator iterViews = sfm_data.views.begin();
    Views::const_iterator iterViewsEnd = sfm_data.views.end();
    if(rangeStart != -1)
    {
      if(rangeStart < 0 || rangeStart > sfm_data.views.size())
      {
        std::cerr << "Bad specific index" << std::endl;
        return EXIT_FAILURE;
      }
      if(rangeSize < 0)
      {
        std::cerr << "Bad range size. " << std::endl;
        return EXIT_FAILURE;
      }
      if(rangeStart + rangeSize > sfm_data.views.size())
        rangeSize = sfm_data.views.size() - rangeStart;

      std::advance(iterViews, rangeStart);
      iterViewsEnd = iterViews;
      std::advance(iterViewsEnd, rangeSize);
    }

    Image<unsigned char> imageGray;
    for(;
      iterViews != iterViewsEnd;
      ++iterViews, ++my_progress_bar)
    {

      const View * view = iterViews->second.get();
      const std::string sView_filename = stlplus::create_filespec(sfm_data.s_root_path,
        view->s_Img_path);
      const std::string sFeat = stlplus::create_filespec(sOutDir,
        stlplus::basename_part(std::to_string(view->id_view)), "feat");
      const std::string sDesc = stlplus::create_filespec(sOutDir,
        stlplus::basename_part(std::to_string(view->id_view)), "desc");

      //If features or descriptors file are missing, compute them
      if (bForce || !stlplus::file_exists(sFeat) || !stlplus::file_exists(sDesc))
      {
        if (!ReadImage(sView_filename.c_str(), &imageGray))
          continue;

        // Compute features and descriptors and export them to files
        std::cout << "Extracting features from image " << view->id_view << std::endl;
        std::unique_ptr<Regions> regions;
        image_describer->Describe(imageGray, regions);
        image_describer->Save(regions.get(), sFeat, sDesc);
      }
    }
    std::cout << "Task done in (s): " << timer.elapsed() << std::endl;
  }
  return EXIT_SUCCESS;
}
