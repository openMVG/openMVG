
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include "openMVG/image/image.hpp"
#include "openMVG/sfm/sfm.hpp"

/// Feature/Regions & Image describer interfaces
#include "openMVG/features/features.hpp"
#include "nonFree/sift/SIFT_describer.hpp"
#include <cereal/archives/json.hpp>
#include "openMVG/system/timer.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "third_party/progress/progress.hpp"

#include <cstdlib>
#include <fstream>

#ifdef OPENMVG_USE_OPENMP
#include <omp.h>
#endif

using namespace openMVG;
using namespace openMVG::image;
using namespace openMVG::features;
using namespace openMVG::sfm;
using namespace std;

features::EDESCRIBER_PRESET stringToEnum(const std::string & sPreset)
{
  features::EDESCRIBER_PRESET preset;
  if(sPreset == "NORMAL")
    preset = features::NORMAL_PRESET;
  else
  if (sPreset == "HIGH")
    preset = features::HIGH_PRESET;
  else
  if (sPreset == "ULTRA")
    preset = features::ULTRA_PRESET;
  else
    preset = features::EDESCRIBER_PRESET(-1);
  return preset;
}

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
#ifdef OPENMVG_USE_OPENMP
  int iNumThreads = 0;
#endif

  // required
  cmd.add( make_option('i', sSfM_Data_Filename, "input_file") );
  cmd.add( make_option('o', sOutDir, "outdir") );
  // Optional
  cmd.add( make_option('m', sImage_Describer_Method, "describerMethod") );
  cmd.add( make_option('u', bUpRight, "upright") );
  cmd.add( make_option('f', bForce, "force") );
  cmd.add( make_option('p', sFeaturePreset, "describerPreset") );

#ifdef OPENMVG_USE_OPENMP
  cmd.add( make_option('n', iNumThreads, "numThreads") );
#endif

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch(const std::string& s) {
      std::cerr << "Usage: " << argv[0] << '\n'
      << "[-i|--input_file] a SfM_Data file \n"
      << "[-o|--outdir path] \n"
      << "\n[Optional]\n"
      << "[-f|--force] Force to recompute data\n"
      << "[-m|--describerMethod]\n"
      << "  (method to use to describe an image):\n"
      << "   SIFT (default),\n"
      << "   AKAZE_FLOAT: AKAZE with floating point descriptors,\n"
      << "   AKAZE_MLDB:  AKAZE with binary descriptors\n"
      << "[-u|--upright] Use Upright feature 0 or 1\n"
      << "[-p|--describerPreset]\n"
      << "  (used to control the Image_describer configuration):\n"
      << "   NORMAL (default),\n"
      << "   HIGH,\n"
      << "   ULTRA: !!Can take long time!!\n"
#ifdef OPENMVG_USE_OPENMP
      << "[-n|--numThreads] number of parallel computations\n"
#endif
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
#ifdef OPENMVG_USE_OPENMP
            << "--numThreads " << iNumThreads << std::endl
#endif
            << std::endl;


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
  if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(VIEWS|INTRINSICS))) {
    std::cerr << std::endl
      << "The input file \""<< sSfM_Data_Filename << "\" cannot be read" << std::endl;
    return false;
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
    // Create the desired Image_describer method.
    // Don't use a factory, perform direct allocation
    if (sImage_Describer_Method == "SIFT")
    {
      image_describer.reset(new SIFT_Image_describer
        (SIFT_Image_describer::Params(), !bUpRight));
    }
    else
    if (sImage_Describer_Method == "AKAZE_FLOAT")
    {
      image_describer.reset(new AKAZE_Image_describer
        (AKAZE_Image_describer::Params(AKAZE::Params(), AKAZE_MSURF), !bUpRight));
    }
    else
    if (sImage_Describer_Method == "AKAZE_MLDB")
    {
      image_describer.reset(new AKAZE_Image_describer
        (AKAZE_Image_describer::Params(AKAZE::Params(), AKAZE_MLDB), !bUpRight));
    }
    if (!image_describer)
    {
      std::cerr << "Cannot create the designed Image_describer:"
        << sImage_Describer_Method << "." << std::endl;
      return EXIT_FAILURE;
    }
    else
    {
      if (!sFeaturePreset.empty())
      if (!image_describer->Set_configuration_preset(stringToEnum(sFeaturePreset)))
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
  // - if regions file exists continue,
  // - if no file, compute features
  {
    system::Timer timer;
    Image<unsigned char> imageGray, globalMask, imageMask;

    const std::string sGlobalMask_filename = stlplus::create_filespec(sOutDir, "mask.png");
    if(stlplus::file_exists(sGlobalMask_filename))
      ReadImage(sGlobalMask_filename.c_str(), &globalMask);

    C_Progress_display my_progress_bar( sfm_data.GetViews().size(),
      std::cout, "\n- EXTRACT FEATURES -\n" );

    #ifdef OPENMVG_USE_OPENMP
    const unsigned int nb_max_thread = omp_get_max_threads();
    #endif

#ifdef OPENMVG_USE_OPENMP
    omp_set_num_threads(iNumThreads);
    #pragma omp parallel for schedule(dynamic) if(iNumThreads > 0) private(imageMask)
#endif
    for(int i = 0; i < sfm_data.views.size(); ++i)
    {
#ifdef OPENMVG_USE_OPENMP
      if(iNumThreads == 0) omp_set_num_threads(nb_max_thread);
#endif
      Views::const_iterator iterViews = sfm_data.views.begin();
      std::advance(iterViews, i);
      const View * view = iterViews->second.get();
      const std::string
        sView_filename = stlplus::create_filespec(sfm_data.s_root_path, view->s_Img_path),
        sFeat = stlplus::create_filespec(sOutDir, stlplus::basename_part(sView_filename), "feat"),
        sDesc = stlplus::create_filespec(sOutDir, stlplus::basename_part(sView_filename), "desc");

      //If features or descriptors file are missing, compute them
      if (bForce || !stlplus::file_exists(sFeat) || !stlplus::file_exists(sDesc))
      {
        if (!ReadImage(sView_filename.c_str(), &imageGray))
          continue;

        Image<unsigned char> * mask = nullptr; // The mask is null by default

        const std::string sImageMask_filename =
          stlplus::create_filespec(sfm_data.s_root_path,
            stlplus::basename_part(sView_filename) + "_mask", "png");

        if(stlplus::file_exists(sImageMask_filename))
          ReadImage(sImageMask_filename.c_str(), &imageMask);

        // The mask point to the globalMask, if a valid one exists for the current image
        if(globalMask.Width() == imageGray.Width() && globalMask.Height() == imageGray.Height())
          mask = &globalMask;
        // The mask point to the imageMask (individual mask) if a valid one exists for the current image
        if(imageMask.Width() == imageGray.Width() && imageMask.Height() == imageGray.Height())
          mask = &imageMask;

        Image<unsigned char> imageGray;
        if (ReadImage(sView_filename.c_str(), &imageGray))
        {
          // Compute features and descriptors and export them to files
          std::unique_ptr<Regions> regions;
          image_describer->Describe(imageGray, regions, mask);
          image_describer->Save(regions.get(), sFeat, sDesc);
        }
      }
#ifdef OPENMVG_USE_OPENMP
      #pragma omp critical
#endif
      ++my_progress_bar;
    }
    std::cout << "Task done in (s): " << timer.elapsed() << std::endl;
  }
  return EXIT_SUCCESS;
}
