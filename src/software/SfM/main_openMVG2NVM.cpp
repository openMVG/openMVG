// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2016 Romuald Perrot, Pierre Moulon

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/cameras/Camera_Pinhole.hpp"
#include "openMVG/cameras/Camera_undistort_image.hpp"
#include "openMVG/features/feature.hpp"
#include "openMVG/image/image_io.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/system/loggerprogress.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <iomanip>
#include <cstdlib>

#ifdef OPENMVG_USE_OPENMP
#include <omp.h>
#endif

using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::geometry;
using namespace openMVG::image;
using namespace openMVG::sfm;
using namespace openMVG::features;

/**
* @brief Create a N-View file and export images
* @param sfm_data Structure from Motion file
* @param sOutDirectory Output directory
* @param filename Name of the file to create
*/
bool CreateNVMFile( const SfM_Data & sfm_data ,
                    const std::string & sOutDirectory ,
                    const std::string & filename )
{
  const std::string sOutViewsDirectory = stlplus::folder_append_separator( sOutDirectory ) + "views";
  if ( !stlplus::folder_exists( sOutViewsDirectory ) )
  {
    OPENMVG_LOG_INFO << "[Creating directory:  " << sOutViewsDirectory << "]";
    stlplus::folder_create( sOutViewsDirectory );
  }

  // Header
  std::ofstream file( filename );

  if ( ! file )
  {
    OPENMVG_LOG_ERROR << "Cannot write file" << filename;
    return false;
  }
  file << "NVM_V3" << "\n";

  // we reindex the poses to ensure a contiguous pose list.
  Hash_Map<IndexT, IndexT> map_viewIdToContiguous;
  int nb_cam = 0;
  for (Views::const_iterator iter = sfm_data.GetViews().begin();
       iter != sfm_data.GetViews().end(); ++iter )
  {
    const View * view = iter->second.get();
    if ( !sfm_data.IsPoseAndIntrinsicDefined( view ) )
    {
      continue;
    }

    nb_cam ++;
    map_viewIdToContiguous.insert( std::make_pair( view->id_view, map_viewIdToContiguous.size() ) );
  }

  // Number of cameras
  // For each camera : File_name Focal Qw Qx Qy Qz Cx Cy Cz D0 0

  file << nb_cam << "\n";

  // Export undistorted images
  {
    system::LoggerProgress my_progress_bar( sfm_data.GetViews().size(), "- EXPORT UNDISTORTED IMAGES -" );
    Image<RGBColor> image, image_ud;
  #ifdef OPENMVG_USE_OPENMP
      #pragma omp parallel for schedule(dynamic) private(image, image_ud)
  #endif
    for (int i = 0; i < static_cast<int>(sfm_data.views.size()); ++i)
    {
      Views::const_iterator iterViews = sfm_data.views.begin();
      std::advance(iterViews, i);
      const View * view = iterViews->second.get();

      if ( !sfm_data.IsPoseAndIntrinsicDefined( view ) )
      {
        continue;
      }

      std::ostringstream padding;
      padding << std::setw( 4 ) << std::setfill( '0' ) << view->id_view;
      const std::string sAbsoluteOutputDir = "view_" + padding.str();
      const std::string sFullOutputDir =
        stlplus::folder_append_separator(sOutViewsDirectory) +
        stlplus::folder_append_separator(sAbsoluteOutputDir);

      const std::string dstImage = stlplus::create_filespec( sFullOutputDir, "undistorted", "png");
      const std::string srcImage = stlplus::create_filespec( sfm_data.s_root_path, view->s_Img_path );

      #ifdef OPENMVG_USE_OPENMP
        #pragma omp critical
      #endif
      // Create output dir if not present
      if ( !stlplus::folder_exists( sFullOutputDir ) )
      {
        stlplus::folder_create( sFullOutputDir );
      }

      Intrinsics::const_iterator iterIntrinsic = sfm_data.GetIntrinsics().find( view->id_intrinsic );
      const IntrinsicBase * cam = iterIntrinsic->second.get();

      // Remove distortion
      if ( cam->have_disto() )
      {
        // Undistort and save the image
        ReadImage( srcImage.c_str(), &image );
        UndistortImage( image, cam, image_ud, BLACK );
        WriteImage( dstImage.c_str(), image_ud );
      }
      else // (no distortion)
      {
        // If extensions match, copy the PNG image
        if ( stlplus::extension_part( srcImage ) == "PNG" ||
             stlplus::extension_part( srcImage ) == "png" )
        {
          stlplus::file_copy( srcImage, dstImage );
        }
        else
        {
          ReadImage( srcImage.c_str(), &image );
          WriteImage( dstImage.c_str(), image );
        }
      }
      ++my_progress_bar;
    }
  }

  // Export camera parameters
  {
    system::LoggerProgress my_progress_bar( sfm_data.GetViews().size(), "- EXPORT CAMERA PARAMETERS -" );
    for (Views::const_iterator iter = sfm_data.GetViews().begin();
         iter != sfm_data.GetViews().end(); ++iter, ++my_progress_bar)
    {
      const View * view = iter->second.get();

      if ( !sfm_data.IsPoseAndIntrinsicDefined( view ) )
      {
        continue;
      }

      std::ostringstream padding;
      padding << std::setw( 4 ) << std::setfill( '0' ) << view->id_view;
      const std::string sAbsoluteOutputDir = stlplus::folder_append_separator("views") + "view_" + padding.str();
      const std::string dstImage = stlplus::create_filespec( sAbsoluteOutputDir, "undistorted", "png");

      Intrinsics::const_iterator iterIntrinsic = sfm_data.GetIntrinsics().find( view->id_intrinsic );
      const IntrinsicBase * cam = iterIntrinsic->second.get();
      const Pinhole_Intrinsic * pinhole_cam = static_cast<const Pinhole_Intrinsic *>( cam );
      const double flen = pinhole_cam->focal();
      const Pose3 pose = sfm_data.GetPoseOrDie( view );
      const Mat3 rotation = pose.rotation();
      const Vec3 center = pose.center();

      const double Cx = center[0];
      const double Cy = center[1];
      const double Cz = center[2];
      Eigen::Quaterniond q( rotation );
      const double Qx = q.x();
      const double Qy = q.y();
      const double Qz = q.z();
      const double Qw = q.w();
      const double d0 = 0.0;

      file << dstImage << " "
         << flen << " "

         << Qw << " "
         << Qx << " "
         << Qy << " "
         << Qz << " "

         << Cx << " "
         << Cy << " "
         << Cz << " "
         << d0 << " "
         << 0 << "\n";
    }
  }

  // Now exports points
  // Number of points
  // For each points : X Y Z R G B Nm [ measurements ]
  // mesurements : Img_idx Feat_idx X Y
  const Landmarks & landmarks = sfm_data.GetLandmarks();
  const size_t featureCount = landmarks.size();
  file << featureCount << "\n";
  system::LoggerProgress my_progress_bar( featureCount, "- EXPORT LANDMARKS DATA -" );
  for ( Landmarks::const_iterator iterLandmarks = landmarks.begin();
        iterLandmarks != landmarks.end(); ++iterLandmarks, ++my_progress_bar )
  {
    const Vec3 exportPoint = iterLandmarks->second.X;
    file << exportPoint.x() << " " << exportPoint.y() << " " << exportPoint.z() << " ";
    file << 250 << " " << 100 << " " << 150 << " ";  // Write arbitrary RGB color

    // Tally set of feature observations
    const Observations & obs = iterLandmarks->second.obs;
    const size_t featureCount = std::distance( obs.begin(), obs.end() );
    file << featureCount;

    for ( Observations::const_iterator itObs = obs.begin(); itObs != obs.end(); ++itObs )
    {
      const IndexT viewId = map_viewIdToContiguous.at(itObs->first);
      const IndexT featId = itObs->second.id_feat;
      const Observation & ob = itObs->second;

      file << " " << viewId << " " << featId << " " << ob.x( 0 ) << " " << ob.x( 1 ) << " ";
    }
    file << "\n";
  }
  // EOF indicator
  file << "0";

  return true;
}

/**
* @brief Main function used to export a NVM file
* @param sfm_data Structure from Motion file to export
* @param sOutDirectory Output directory
*/
bool exportToNVM( const SfM_Data & sfm_data , const std::string & sOutDirectory  )
{
  // Create output directory
  bool bOk = false;
  if ( !stlplus::is_folder( sOutDirectory ) )
  {
    OPENMVG_LOG_INFO << "[Creating directory:  " << sOutDirectory << "]";
    stlplus::folder_create( sOutDirectory );
    bOk = stlplus::is_folder( sOutDirectory );
  }
  else
  {
    bOk = true;
  }

  if ( !bOk )
  {
    OPENMVG_LOG_ERROR << "Cannot access one of the desired output directories";
    return false;
  }
  const std::string sFilename = stlplus::create_filespec( sOutDirectory , "scene.nvm" );
  if ( ! CreateNVMFile( sfm_data , sOutDirectory , sFilename ) )
  {
    OPENMVG_LOG_ERROR << "There was an error exporting project";
    return false;
  }
  return true;

}

int main( int argc , char ** argv )
{
  CmdLine cmd;
  std::string sSfM_Data_Filename;
  std::string sOutDir = "";
#ifdef OPENMVG_USE_OPENMP
  int iNumThreads = 1;
#endif

  cmd.add( make_option( 'i', sSfM_Data_Filename, "sfmdata" ) );
  cmd.add( make_option( 'o', sOutDir, "outdir" ) );
#ifdef OPENMVG_USE_OPENMP
  cmd.add( make_option('n', iNumThreads, "numThreads") );
#endif
  OPENMVG_LOG_INFO << "Note:  this program writes output in NVM file format.";

  try
  {
    if ( argc == 1 )
    {
      throw std::string( "Invalid command line parameter." );
    }
    cmd.process( argc, argv );
  }
  catch ( const std::string& s )
  {
    OPENMVG_LOG_INFO << "Usage: " << argv[0] << '\n'
      << "[-i|--sfmdata] filename, the SfM_Data file to convert\n"
      << "[-o|--outdir] path where the scene.nvm will be saved\n"
#ifdef OPENMVG_USE_OPENMP
      << "[-n|--numThreads] number of thread(s)"
#endif
      ;

    OPENMVG_LOG_ERROR << s;
    return EXIT_FAILURE;
  }

  if ( !stlplus::folder_exists( sOutDir ) )
  {
    stlplus::folder_create( sOutDir );
  }

#ifdef OPENMVG_USE_OPENMP
    const unsigned int nb_max_thread = omp_get_max_threads();

    if (iNumThreads > 0) {
        omp_set_num_threads(iNumThreads);
    } else {
        omp_set_num_threads(nb_max_thread);
    }
#endif

  // Read the input SfM scene
  SfM_Data sfm_data;
  if ( !Load( sfm_data, sSfM_Data_Filename, ESfM_Data( ALL ) ) )
  {
    OPENMVG_LOG_ERROR << "The input SfM_Data file \"" << sSfM_Data_Filename << "\" cannot be read.";
    return EXIT_FAILURE;
  }

  if ( ! exportToNVM( sfm_data , sOutDir ) )
  {
    OPENMVG_LOG_ERROR << "There was an error during exporting the NVM the file";
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
