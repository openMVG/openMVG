// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/cameras/Camera_undistort_image.hpp"
#include "openMVG/image/image_io.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/system/timer.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/progress/progress_display.hpp"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <cstdlib>
#include <string>

#ifdef OPENMVG_USE_OPENMP
#include <omp.h>
#endif

using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::geometry;
using namespace openMVG::image;
using namespace openMVG::sfm;

#ifdef OPENMVG_USE_OPENMP
#include <omp.h>
#endif

int main(int argc, char *argv[]) {

  CmdLine cmd;
  std::string sSfM_Data_Filename;
  std::string sOutDir = "";
  bool bExportOnlyReconstructedViews = false;
#ifdef OPENMVG_USE_OPENMP
  int iNumThreads = 0;
#endif

  cmd.add( make_option('i', sSfM_Data_Filename, "sfmdata") );
  cmd.add( make_option('o', sOutDir, "outdir") );
  cmd.add( make_option('r', bExportOnlyReconstructedViews, "exportOnlyReconstructed") );

#ifdef OPENMVG_USE_OPENMP
  cmd.add( make_option('n', iNumThreads, "numThreads") );
#endif

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch (const std::string& s) {
    std::cerr
      << "Export undistorted images related to a sfm_data file.\n"
      << "Usage: " << argv[0] << '\n'
      << "[-i|--sfmdata] filename, the SfM_Data file to convert\n"
      << "[-o|--outdir] path\n"
      << "[-r|--exportOnlyReconstructed] boolean 1/0 (default = 0)\n"
#ifdef OPENMVG_USE_OPENMP
      << "[-n|--numThreads] number of thread(s)\n"
#endif
      << std::endl;

      std::cerr << s << std::endl;
      return EXIT_FAILURE;
  }

  // Create output dir
  if (!stlplus::folder_exists(sOutDir))
    stlplus::folder_create( sOutDir );

  SfM_Data sfm_data;
  if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(VIEWS|INTRINSICS|EXTRINSICS))) {
    std::cerr << std::endl
      << "The input SfM_Data file \""<< sSfM_Data_Filename << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }

  bool bOk = true;
  {
    system::Timer timer;
    // Export views as undistorted images (those with valid Intrinsics)
    Image<RGBColor> image, image_ud;
    Image<uint8_t> image_gray, image_gray_ud;
    C_Progress_display my_progress_bar( sfm_data.GetViews().size(), std::cout, "\n- EXTRACT UNDISTORTED IMAGES -\n" );

    #ifdef OPENMVG_USE_OPENMP
    const unsigned int nb_max_thread = omp_get_max_threads();
    #endif

#ifdef OPENMVG_USE_OPENMP
    omp_set_num_threads(iNumThreads);
    #pragma omp parallel for schedule(dynamic) if (iNumThreads > 0) private(image, image_ud, image_gray, image_gray_ud)
#endif
    for (int i = 0; i < static_cast<int>(sfm_data.views.size()); ++i)
    {
#ifdef OPENMVG_USE_OPENMP
      if (iNumThreads == 0) omp_set_num_threads(nb_max_thread);
#endif
      Views::const_iterator iterViews = sfm_data.views.begin();
      std::advance(iterViews, i);

      const View * view = iterViews->second.get();
      // Check if the view is in reconstruction
      if (bExportOnlyReconstructedViews && !sfm_data.IsPoseAndIntrinsicDefined(view))
        continue;

      const bool bIntrinsicDefined = view->id_intrinsic != UndefinedIndexT &&
        sfm_data.GetIntrinsics().find(view->id_intrinsic) != sfm_data.GetIntrinsics().end();
      if (!bIntrinsicDefined)
        continue;

      Intrinsics::const_iterator iterIntrinsic = sfm_data.GetIntrinsics().find(view->id_intrinsic);

      const std::string srcImage = stlplus::create_filespec(sfm_data.s_root_path, view->s_Img_path);
      const std::string dstImage = stlplus::create_filespec(
        sOutDir, stlplus::filename_part(srcImage));

      const IntrinsicBase * cam = iterIntrinsic->second.get();
      if (cam->have_disto())
      {
        // undistort the image and save it
        if (ReadImage( srcImage.c_str(), &image))
        {
          UndistortImage(image, cam, image_ud, BLACK);
          const bool bRes = WriteImage(dstImage.c_str(), image_ud);
#ifdef OPENMVG_USE_OPENMP
          #pragma omp critical
#endif
          bOk &= bRes;
        }
        else // If RGBColor reading fails, we try to read a gray image
        if (ReadImage( srcImage.c_str(), &image_gray))
        {
          UndistortImage(image_gray, cam, image_gray_ud, BLACK);
          const bool bRes = WriteImage(dstImage.c_str(), image_gray_ud);
#ifdef OPENMVG_USE_OPENMP
          #pragma omp critical
#endif
          bOk &= bRes;
        }
      }
      else // (no distortion)
      {
        // copy the image since there is no distortion
        stlplus::file_copy(srcImage, dstImage);
      }
      ++my_progress_bar;
    }
    std::cout << "Task done in (s): " << timer.elapsed() << std::endl;
  }

  // Exit program
  if (bOk)
  {
    return EXIT_SUCCESS;
  }
  return EXIT_FAILURE;
}
