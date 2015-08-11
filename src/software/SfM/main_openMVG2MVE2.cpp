/* v.0.15 7 August 2015
 * Kevin CAIN, www.insightdigital.org
 * Adapted from the openMVG libraries,
 * Copyright (c) 2012-2015 Pierre MOULON.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "openMVG/sfm/sfm.hpp"
#include "openMVG/image/image.hpp"

using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::geometry;
using namespace openMVG::image;
using namespace openMVG::sfm;
using namespace openMVG::features;

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/progress/progress.hpp"
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <iterator>
#include <iomanip>

/* Notes:
 * - An MVE2 scene appears to duplicate camera rot matrix and trans vector per-view data in 'meta.ini'
 *   within the first section of 'synth_0.out'.
 * - We do not save the original, instead we rely on the undistorted image from openMVG.
 * - We do not output thumbnails or EXIF blobs, as these appear only to be used only for the GUI UMVE.
 * - To avoid encoding loss, openMVG images should be written as .PNG if undistorted images are *not* computed.
 * - In OpenMVG, some views may have some missing poses; MVE does *not* require a contiguous camera index.
 *
 *  For information on the target for this conversion, please see the MVE (v2) File format:
 *  https://github.com/simonfuhrmann/mve/wiki/MVE-File-Format
 */

bool exportToMVE2Format(
  const SfM_Data & sfm_data,
  const std::string & sOutDirectory // Output MVE2 files directory
  )
{
  bool bOk = true;
  // Create basis directory structure
  if (!stlplus::is_folder(sOutDirectory))
  {
    cout << "\033[1;31mCreating directory:  " << sOutDirectory << "\033[0m\n";
    stlplus::folder_create(sOutDirectory);
    bOk = stlplus::is_folder(sOutDirectory);
  }

  if (!bOk)
  {
    std::cerr << "Cannot access one of the desired output directories" << std::endl;
	  return false;
  }
  
  // Export the SfM_Data scene to the MVE2 format
  {
    // Create 'views' subdirectory
    const string sOutViewsDirectory = stlplus::folder_append_separator(sOutDirectory) + "views";
    if (!stlplus::folder_exists(sOutViewsDirectory))
    {
      cout << "\033[1;31mCreating directory:  " << sOutViewsDirectory << "\033[0m\n";
      stlplus::folder_create(sOutViewsDirectory);
    }

    // Prepare to write bundle file
    // Get cameras and features from OpenMVG
    const size_t cameraCount = std::distance(sfm_data.GetViews().begin(), sfm_data.GetViews().end());
    // Tally global set of feature landmarks
    const Landmarks & landmarks = sfm_data.GetLandmarks();
    const size_t featureCount = std::distance(landmarks.begin(), landmarks.end());
    const std::string filename = "synth_0.out";
    std::cout << "Writing bundle (" << cameraCount << " cameras, "
        << featureCount << " features): to " << filename << "...\n";
    std::ofstream out(stlplus::folder_append_separator(sOutDirectory) + filename);
    out << "drews 1.0\n";  // MVE expects this header
    out << cameraCount << " " << featureCount << "\n";

    // Export (calibrated) views as undistorted images
    C_Progress_display my_progress_bar( sfm_data.GetViews().size()*1);
    std::pair<int,int> w_h_image_size;
    Image<RGBColor> image, image_ud;
    std::string sOutViewIteratorDirectory;
    for(Views::const_iterator iter = sfm_data.GetViews().begin();
      iter != sfm_data.GetViews().end(); ++iter, ++my_progress_bar)
    {
      const View * view = iter->second.get();
      if (!sfm_data.IsPoseAndIntrinsicDefined(view))
        continue;

      // Create current view subdirectory 'view_xxxx.mve'
      std::ostringstream padding;
      padding << std::setw(4) << std::setfill('0') << view->id_view;
      sOutViewIteratorDirectory = stlplus::folder_append_separator(sOutViewsDirectory) + "view_" + padding.str() + ".mve";
      if (!stlplus::folder_exists(sOutViewIteratorDirectory))
      {
        stlplus::folder_create(sOutViewIteratorDirectory);
      }
      Intrinsics::const_iterator iterIntrinsic = sfm_data.GetIntrinsics().find(view->id_intrinsic);

      // We have a valid view with a corresponding camera & pose
      const std::string srcImage = stlplus::create_filespec(sfm_data.s_root_path, view->s_Img_path);
      std::string dstImage =
        stlplus::create_filespec(stlplus::folder_append_separator(sOutViewIteratorDirectory), "undistorted","png");

      const IntrinsicBase * cam = iterIntrinsic->second.get();
      if (cam->have_disto())
      {
        // Undistort and save the image
        ReadImage(srcImage.c_str(), &image);
        UndistortImage(image, cam, image_ud, BLACK);
        WriteImage(dstImage.c_str(), image_ud);
      }
      else // (no distortion)
      {
        // If extensions match, copy the PNG image
        if (stlplus::extension_part(srcImage) == "PNG" ||
          stlplus::extension_part(srcImage) == "png")
        {
          stlplus::file_copy(srcImage, dstImage);
        }
        else
        {
          ReadImage( srcImage.c_str(), &image);
          WriteImage( dstImage.c_str(), image);
        }
      }

      // Prepare to write an MVE 'meta.ini' file for the current view
      const Pose3 pose = sfm_data.GetPoseOrDie(view);
      Mat3 K; // intrinsic camera parameter matrix      
      const Pinhole_Intrinsic * pinhole_cam = static_cast<const Pinhole_Intrinsic *>(cam);
      if (pinhole_cam)
      {
        K = pinhole_cam->K();
      }
      else
      {
        const Mat34 P = cam->get_projective_equivalent(pose);
        Mat3 R;
        Vec3 t;
        KRt_From_P(P, &K, &R, &t);
      }

      const Mat3 rotation = pose.rotation();
      const Vec3 translation = pose.translation();
      // Pixel aspect: assuming square pixels
      const float pixelAspect = 1.f;
      // Focal length and principal point are embedded within calibration matrix K:
      // focal_length = K(0,0)
      // principal point = {_K(0,2), _K(1,2)}
      // their values are normalized (0..1)
      const float flen = K(0,0) / static_cast<double>(std::max(cam->w(), cam->h()));
      const float ppX = std::abs(K(0,2)/cam->w());
      const float ppY = std::abs(K(1,2)/cam->h());

      std::ostringstream fileOut;
      fileOut << "#MVE view meta data is stored in INI-file syntax." << fileOut.widen('\n')
        << "#This file is generated, formatting will get lost." << fileOut.widen('\n')
        << fileOut.widen('\n')
        << "[camera]" << fileOut.widen('\n')
        << "focal_length = " << flen << fileOut.widen('\n')
        << "pixel_aspect = " << pixelAspect << fileOut.widen('\n')
        << "principal_point = " << ppX << " " << ppY << fileOut.widen('\n')
        << "rotation = " << rotation(0, 0) << " " << rotation(0, 1) << " " << rotation(0, 2) << " "
        << rotation(1, 0) << " " << rotation(1, 1) << " " << rotation(1, 2) << " "
        << rotation(2, 0) << " " << rotation(2, 1) << " " << rotation(2, 2) << fileOut.widen('\n')
        << "translation = " << translation[0] << " " << translation[1] << " "
        << translation[2] << " " << fileOut.widen('\n')
        << "[view]" << fileOut.widen('\n')
        << "id = " << view->id_view << fileOut.widen('\n')
        << "name = " << srcImage.c_str() << fileOut.widen('\n');

      // To do:  trim any extra separator(s) from openMVG name we receive, e.g.:
      // '/home/insight/openMVG_KevinCain/openMVG_Build/software/SfM/ImageDataset_SceauxCastle/images//100_7100.JPG'
      std::ofstream file(
  	    stlplus::create_filespec(stlplus::folder_append_separator(sOutViewIteratorDirectory),
	      "meta","ini").c_str());
      file << fileOut.str();
      file.close();

      // Save a thumbnail image "thumbnail.png", 50x50 pixels
      // For now, we ignore thumbnails under the assumption that they are used only by UMVE
      // Pierre Moulon suggested we resample as per:  https://github.com/openMVG/openMVG/blob/develop/src/openMVG/image/image_resampling_test.cpp#L24

      // For each camera, write to bundle:  focal length, radial distortion[0-1], rotation matrix[0-8], translation vector[0-2]
      // To do:  add more rigorous camera sanity checks, as per:
      // https://github.com/simonfuhrmann/mve/blob/e3db7bc60ce93fe51702ba77ef480e151f927c23/libs/mve/bundle_io.cc
      if (flen == 0.0f)
      {
        for (int i = 0; i < 5 * 3; ++i)
          out << "0" << (i % 3 == 2 ? "\n" : " ");
        continue;
      }

      out << flen << " " << "0" << " " << "0" << "\n"  // Write '0' distortion values for pre-corrected images
        << rotation(0, 0) << " " << rotation(0, 1) << " " << rotation(0, 2) << "\n"
        << rotation(1, 0) << " " << rotation(1, 1) << " " << rotation(1, 2) << "\n"
        << rotation(2, 0) << " " << rotation(2, 1) << " " << rotation(2, 2) << "\n"
        << translation[0] << " " << translation[1] << " " << translation[2] << "\n";
    }

    // For each feature, write to bundle:  position XYZ[0-3], color RGB[0-2], all ref.view_id & ref.feature_id
    // The following method is adapted from Simon Fuhrmann's MVE project:
    // https://github.com/simonfuhrmann/mve/blob/e3db7bc60ce93fe51702ba77ef480e151f927c23/libs/mve/bundle_io.cc

    for (Landmarks::const_iterator iterLandmarks = landmarks.begin(); iterLandmarks != landmarks.end(); ++iterLandmarks)  {
      const Vec3 exportPoint = iterLandmarks->second.X;
      out << exportPoint.x() << " " << exportPoint.y() << " " << exportPoint.z() << "\n";
      out << 250 << " " << 100 << " " << 150 << "\n";  // Write arbitrary RGB color, see above note

      // Tally set of feature observations
      const Observations & obs = iterLandmarks->second.obs;
      const size_t featureCount = std::distance(obs.begin(), obs.end());
      out << featureCount;

      for (Observations::const_iterator itObs = obs.begin(); itObs != obs.end(); ++itObs)  {
          const IndexT viewId = itObs->first;
          const IndexT featId = itObs->second.id_feat;
          out << " " << viewId << " " << featId << " 0";
      }
      out << "\n";
    }
    out.close();
  }
  return bOk;
}

int main(int argc, char *argv[]) {

  CmdLine cmd;
  std::string sSfM_Data_Filename;
  std::string sOutDir = "";
  cmd.add( make_option('i', sSfM_Data_Filename, "sfmdata") );
  cmd.add( make_option('o', sOutDir, "outdir") );
  std::cout << "Note:  this program writes output in MVE file format.\n";

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch(const std::string& s) {
      std::cerr << "Usage: " << argv[0] << '\n'
      << "[-i|--sfmdata] filename, the SfM_Data file to convert\n"
      << "[-o|--outdir] path\n"
      << std::endl;

      std::cerr << s << std::endl;
      return EXIT_FAILURE;
  }

  // Create output dir
  if (!stlplus::folder_exists(sOutDir))
    stlplus::folder_create(sOutDir);

  // Read the input SfM scene
  SfM_Data sfm_data;
  if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(ALL))) {
    std::cerr << std::endl
      << "The input SfM_Data file \""<< sSfM_Data_Filename << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }

  if (exportToMVE2Format(sfm_data, stlplus::folder_append_separator(sOutDir) + "MVE"))
    return( EXIT_SUCCESS );
  else
    return( EXIT_FAILURE );
}
