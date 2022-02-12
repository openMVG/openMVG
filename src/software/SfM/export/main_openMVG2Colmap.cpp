// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2016 Romuald Perrot, Pierre Moulon

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/cameras/Camera_Common.hpp"
#include "openMVG/cameras/Camera_Pinhole.hpp"
#include "openMVG/cameras/Camera_Pinhole_Radial.hpp"
#include "openMVG/cameras/Camera_Pinhole_Fisheye.hpp"
#include "openMVG/cameras/Camera_Intrinsics.hpp"
#include "openMVG/cameras/Camera_undistort_image.hpp"
#include "openMVG/features/feature.hpp"
#include "openMVG/image/image_io.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_data_colorization.hpp"
#include "openMVG/system/loggerprogress.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <iomanip>
#include <cstdlib>
#include <sstream>

using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::geometry;
using namespace openMVG::image;
using namespace openMVG::sfm;
using namespace openMVG::features;


bool CreateLineCameraFile(  const IndexT camera_id, 
                            std::shared_ptr<openMVG::cameras::IntrinsicBase> intrinsic,
                            std::string & camera_linie,
                            const int& floating_point_precision_digit)

{
  std::stringstream came_line_ss;
  came_line_ss.precision(floating_point_precision_digit);

  EINTRINSIC current_type = intrinsic->getType();
  switch(current_type) 
  {
    case PINHOLE_CAMERA: 
      //OpenMVG's PINHOLE_CAMERA corresponds to Colmap's SIMPLE_PINHOLE
      //Parameters: f, cx, cy  
      {
        std::shared_ptr<openMVG::cameras::Pinhole_Intrinsic> pinhole_intrinsic(
          dynamic_cast<openMVG::cameras::Pinhole_Intrinsic * >(intrinsic->clone()));
        
        came_line_ss << camera_id << " " << 
          "SIMPLE_PINHOLE" << " " <<
          pinhole_intrinsic->w() << " " << 
          pinhole_intrinsic->h() << " " <<
          pinhole_intrinsic->focal() << " " << 
          pinhole_intrinsic->principal_point().x() << " " << 
          pinhole_intrinsic->principal_point().y() << "\n";
      }
      break;
    case PINHOLE_CAMERA_RADIAL1:
      //OpenMVG's PINHOLE_CAMERA_RADIAL1 corresponds to Colmap's SIMPLE_RADIAL
      //Parameters: f, cx, cy, k1   
      {
        std::shared_ptr<openMVG::cameras::Pinhole_Intrinsic_Radial_K1> pinhole_intrinsic_radial(
            dynamic_cast<openMVG::cameras::Pinhole_Intrinsic_Radial_K1 * >(intrinsic->clone()));

        came_line_ss << camera_id << " " << 
          "SIMPLE_RADIAL" << " " <<
          pinhole_intrinsic_radial->w() << " " << 
          pinhole_intrinsic_radial->h() << " " <<
          pinhole_intrinsic_radial->focal() << " " << 
          pinhole_intrinsic_radial->principal_point().x() << " " << 
          pinhole_intrinsic_radial->principal_point().y() << " " << 
          pinhole_intrinsic_radial->getParams().at(3) << "\n";   //k1
      }
      break;      
    case PINHOLE_CAMERA_RADIAL3: 
      //OpenMVG's PINHOLE_CAMERA_RADIAL3 corresponds to Colmap's FULL_OPENCV
      //Parameters: fx, fy, cx, cy, k1, k2, p1, p2, k3, k4, k5, k6
      {
        std::shared_ptr<openMVG::cameras::Pinhole_Intrinsic_Radial_K3> pinhole_intrinsic_radial(
            dynamic_cast<openMVG::cameras::Pinhole_Intrinsic_Radial_K3 * >(intrinsic->clone()));
        
        came_line_ss << camera_id << " " << 
          "FULL_OPENCV" << " " <<
          pinhole_intrinsic_radial->w() << " " << 
          pinhole_intrinsic_radial->h() << " " <<
          pinhole_intrinsic_radial->focal() << " " << 
          pinhole_intrinsic_radial->focal() << " " << 
          pinhole_intrinsic_radial->principal_point().x() << " " << 
          pinhole_intrinsic_radial->principal_point().y() << " " << 
          pinhole_intrinsic_radial->getParams().at(3) << " " << //k1
          pinhole_intrinsic_radial->getParams().at(4) << " " << //k2
          0.0 << " " << //p1
          0.0 << " " << //p2
          pinhole_intrinsic_radial->getParams().at(5) << " " << //k3
          0.0 << " " << // k4
          0.0 << " " << // k5
          0.0 << "\n";  // k6
      }
      break; 
    case PINHOLE_CAMERA_BROWN: 
      OPENMVG_LOG_ERROR << "PINHOLE_CAMERA_BROWN is not supported. Aborting ...";
      return false;
      break;      
    case PINHOLE_CAMERA_FISHEYE: 
      //OpenMVG's PINHOLE_CAMERA_FISHEYE corresponds to Colmap's OPENCV_FISHEYE
      //Parameters: fx, fy, cx, cy, k1, k2, k3, k4
      {
        std::shared_ptr<openMVG::cameras::Pinhole_Intrinsic_Fisheye> pinhole_intrinsic_fisheye(
            dynamic_cast<openMVG::cameras::Pinhole_Intrinsic_Fisheye * >(intrinsic->clone()));
        
        came_line_ss << camera_id << " " << 
          "OPENCV_FISHEYE" << " " <<
          pinhole_intrinsic_fisheye->w() << " " << 
          pinhole_intrinsic_fisheye->h() << " " <<
          pinhole_intrinsic_fisheye->focal() << " " << 
          pinhole_intrinsic_fisheye->focal() << " " << 
          pinhole_intrinsic_fisheye->principal_point().x() << " " << 
          pinhole_intrinsic_fisheye->principal_point().y() << " " << 
          pinhole_intrinsic_fisheye->getParams().at(3) << " " << //k1
          pinhole_intrinsic_fisheye->getParams().at(4) << " " << //k2
          pinhole_intrinsic_fisheye->getParams().at(5) << " " << //k3
          pinhole_intrinsic_fisheye->getParams().at(6) << " " << "\n"; //k4
      }
      break; 
    default: OPENMVG_LOG_ERROR << "Camera Type " << current_type << " is not supported. Aborting ...";
    return false;
  }
  camera_linie = came_line_ss.str();
  return true;
}

bool CreateCameraFile( const SfM_Data & sfm_data, 
                      const std::string & sCamerasFilename,
                      const int& floating_point_precision_digit)
{
   /* cameras.txt
      # Camera list with one line of data per camera:
      #   CAMERA_ID, MODEL, WIDTH, HEIGHT, PARAMS[]
      # Number of cameras: X
  */ 
  std::ofstream camera_file( sCamerasFilename );
  camera_file.precision(floating_point_precision_digit);

  if ( ! camera_file )
  {
    OPENMVG_LOG_ERROR << "Cannot write file" << sCamerasFilename;
    return false;
  }
  camera_file << "# Camera list with one line of data per camera:\n";
  camera_file << "#   CAMERA_ID, MODEL, WIDTH, HEIGHT, PARAMS[]\n";
  camera_file << "# Number of cameras: X\n";

  std::vector<std::string> camera_lines;
  system::LoggerProgress my_progress_bar( sfm_data.GetIntrinsics().size(), "- CREATE CAMERA FILE -" );
  for (Intrinsics::const_iterator iter = sfm_data.GetIntrinsics().begin();
    iter != sfm_data.GetIntrinsics().end(); ++iter, ++my_progress_bar)
  {
    const IndexT camera_id = iter->first;
    std::shared_ptr<openMVG::cameras::IntrinsicBase> intrinsic = iter->second;
    std::string camera_line;
    if (CreateLineCameraFile(camera_id, intrinsic, camera_line, floating_point_precision_digit)) 
    {
      camera_lines.push_back(camera_line);
    } 
    else 
    {
      return false;
    }
  }
  for (auto const& camera_line: camera_lines)
  {
    camera_file << camera_line << "\n";
  }
  return true;
}

bool CreateImageFile( const SfM_Data & sfm_data,
                      const std::string & sImagesFilename,
                      const int& floating_point_precision_digit)

{
 /* images.txt
      # Image list with two lines of data per image:
      #   IMAGE_ID, QW, QX, QY, QZ, TX, TY, TZ, CAMERA_ID, NAME
      #   POINTS2D[] as (X, Y, POINT3D_ID)
      # Number of images: X, mean observations per image: Y
  */

  // Header
  std::ofstream images_file( sImagesFilename );

  if ( ! images_file )
  {
    OPENMVG_LOG_ERROR << "Cannot write file" << sImagesFilename;
    return false;
  }
  images_file << "# Image list with two lines of data per image:\n";
  images_file << "#   IMAGE_ID, QW, QX, QY, QZ, TX, TY, TZ, CAMERA_ID, NAME\n";
  images_file << "#   POINTS2D[] as (X, Y, POINT3D_ID)\n";
  images_file << "# Number of images: X, mean observations per image: Y\n";

  std::map< IndexT, std::vector< std::tuple<double, double, IndexT> > > viewIdToPoints2D;
  const Landmarks & landmarks = sfm_data.GetLandmarks();
  {
    for ( Landmarks::const_iterator iterLandmarks = landmarks.begin();
          iterLandmarks != landmarks.end(); ++iterLandmarks)
    {
      const IndexT point3d_id = iterLandmarks->first;

      // Tally set of feature observations
      const Observations & obs = iterLandmarks->second.obs;
      for ( Observations::const_iterator itObs = obs.begin(); itObs != obs.end(); ++itObs )
      {
        const IndexT currentViewId = itObs->first;
        const Observation & ob = itObs->second;
        viewIdToPoints2D[currentViewId].push_back(std::make_tuple(ob.x( 0 ), ob.x( 1 ), point3d_id));
      }
    }
  }

  {
    system::LoggerProgress my_progress_bar( sfm_data.GetViews().size(), "- CREATE IMAGE FILE -" );

    for (Views::const_iterator iter = sfm_data.GetViews().begin();
         iter != sfm_data.GetViews().end(); ++iter, ++my_progress_bar)
    {
      const View * view = iter->second.get();

      if ( !sfm_data.IsPoseAndIntrinsicDefined( view ) )
      {
        continue;
      }

      const Pose3 pose = sfm_data.GetPoseOrDie( view );
      const Mat3 rotation = pose.rotation();
      const Vec3 translation = pose.translation();

      const double Tx = translation[0];
      const double Ty = translation[1];
      const double Tz = translation[2];
      Eigen::Quaterniond q( rotation );
      const double Qx = q.x();
      const double Qy = q.y();
      const double Qz = q.z();
      const double Qw = q.w();

      const IndexT image_id = view->id_view;
      // Colmap's camera_ids correspond to openMVG's intrinsic ids
      const IndexT camera_id = view->id_intrinsic;                           
      const std::string image_name = view->s_Img_path;

      // first line per image
      //IMAGE_ID, QW, QX, QY, QZ, TX, TY, TZ, CAMERA_ID, NAME
      images_file << image_id << " "
         << Qw << " "
         << Qx << " "
         << Qy << " "
         << Qz << " "
         << Tx << " "
         << Ty << " "
         << Tz << " "
         << camera_id << " "
         << image_name << " "     
         << "\n";

      // second line per image 
      //POINTS2D[] as (X, Y, POINT3D_ID)
      for (auto point2D: viewIdToPoints2D[image_id]) 
      {
        images_file << std::get<0>(point2D) << " " << 
        std::get<1>(point2D) << " " <<
        std::get<2>(point2D) << " ";
      }
      images_file << "\n";
    }
  }
  return true;
}

bool CreatePoint3DFile( const SfM_Data & sfm_data,
                      const std::string & sPoints3DFilename,
                      const int& floating_point_precision_digit)
{
 /* points3D.txt
      # 3D point list with one line of data per point:
      #   POINT3D_ID, X, Y, Z, R, G, B, ERROR, TRACK[] as (IMAGE_ID, POINT2D_IDX)
      # Number of points: X, mean track length: Y
  */

  std::ofstream points3D_file( sPoints3DFilename );
  points3D_file.precision(floating_point_precision_digit);

  if ( ! points3D_file )
  {
    OPENMVG_LOG_ERROR << "Cannot write file" << sPoints3DFilename;
    return false;
  }
  points3D_file << "# 3D point list with one line of data per point:\n";
  points3D_file << "#   POINT3D_ID, X, Y, Z, R, G, B, ERROR, TRACK[] as (IMAGE_ID, POINT2D_IDX)\n";
  points3D_file << "# Number of points: X, mean track length: Y\n";

  const Landmarks & landmarks = sfm_data.GetLandmarks();

  std::vector<Vec3> vec_3dPoints, vec_tracksColor;
  if (!ColorizeTracks(sfm_data, vec_3dPoints, vec_tracksColor)) {
    return false;
  }

  system::LoggerProgress my_progress_bar( landmarks.size(), "- CREATE POINT3D FILE  -" );
  int point_index = 0;
  for ( Landmarks::const_iterator iterLandmarks = landmarks.begin();
        iterLandmarks != landmarks.end(); ++iterLandmarks, ++my_progress_bar )
  {
    const Vec3 exportPoint = iterLandmarks->second.X;
    const IndexT point3d_id = iterLandmarks->first;
    points3D_file << point3d_id << " "
      << exportPoint.x() << " " 
      << exportPoint.y() << " " 
      << exportPoint.z() << " "
    
      << static_cast<int>(vec_tracksColor.at(point_index)(0)) << " " 
      << static_cast<int>(vec_tracksColor.at(point_index)(1)) << " " 
      << static_cast<int>(vec_tracksColor.at(point_index)(2)) << " ";

    ++point_index;

    const double error = 0.0;     // Some error
    points3D_file << error;

    const Observations & obs = iterLandmarks->second.obs;
    for ( Observations::const_iterator itObs = obs.begin(); itObs != obs.end(); ++itObs )
    {
      const IndexT viewId = itObs->first;
      const IndexT featId = itObs->second.id_feat;

      points3D_file << " " 
      << viewId << " " 
      << featId;
    }
    points3D_file << "\n";
  }

  return true;
}

/**
* @brief Convert OpenMVG reconstruction result to Colmap's reconstruction format. 
* @param sfm_data Structure from Motion file
* @param sOutDirectory Output directory
* @param sCamerasFilename File name of the camera file
* @param sImagesFilename File name of the image file
* @param sPoints3DFilename File name of the point3D file
* @param floating_point_precision_digit Decimal precision to be used to format floating-point values for output files
*/
bool CreateColmapFolder( const SfM_Data & sfm_data,
                    const std::string & sOutDirectory,
                    const std::string & sCamerasFilename,
                    const std::string & sImagesFilename, 
                    const std::string & sPoints3DFilename,
                    const int& floating_point_precision_digit)
{
  /* Colmap Output structure:
      cameras.txt
      images.txt
      points3D.txt
  */
  if (!CreateCameraFile(sfm_data, sCamerasFilename, floating_point_precision_digit)) 
  {
    return false;
  }
  if (!CreateImageFile(sfm_data, sImagesFilename, floating_point_precision_digit))
  {
    return false;
  }
  if (! CreatePoint3DFile(sfm_data, sPoints3DFilename, floating_point_precision_digit))
  {
    return false;
  }
  return true;
}

/**
* @brief Main function used to export a Colmap reconstruction folder
* @param sfm_data Structure from Motion file to export
* @param sOutDirectory Output directory
*/
bool exportToColmap( const SfM_Data & sfm_data , const std::string & sOutDirectory , const int& floating_point_precision_digit)
{
  // Create output directory
  bool bOk = false;
  if ( !stlplus::is_folder( sOutDirectory ) )
  {
    OPENMVG_LOG_INFO << "\033[1;31mCreating directory:  " << sOutDirectory << "\033[0m\n";
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
  const std::string sCamerasFilename = stlplus::create_filespec( sOutDirectory , "cameras.txt" );
  const std::string sImagesFilename = stlplus::create_filespec( sOutDirectory , "images.txt" );
  const std::string sPoints3DFilename = stlplus::create_filespec( sOutDirectory , "points3D.txt" );
  if ( ! CreateColmapFolder( sfm_data , sOutDirectory , sCamerasFilename, sImagesFilename, sPoints3DFilename, floating_point_precision_digit) )
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
  int floating_point_precision_digit = 16;

  cmd.add( make_option( 'i', sSfM_Data_Filename, "sfmdata" ) );
  cmd.add( make_option( 'o', sOutDir, "outdir" ) );
  cmd.add(make_option('p', floating_point_precision_digit, "precision"));

  OPENMVG_LOG_INFO << "Note: this program writes output in Colmap file format.\n";

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
    OPENMVG_LOG_INFO 
              << "Usage: " << argv[0] << '\n'
              << "[-i|--sfmdata] filename, the SfM_Data file to convert\n"
              << "[-o|--outdir] path where cameras.txt, images.txt and points3D.txt will be saved"
			  << "\n[Optional]\n"
      		  << "[-p|--precision] sets the decimal precision to be used to format floating-point values (default = 16)";
    OPENMVG_LOG_ERROR << s;
    return EXIT_FAILURE;
  }

  if ( !stlplus::folder_exists( sOutDir ) )
  {
    stlplus::folder_create( sOutDir );
  }

  // Read the input SfM scene
  SfM_Data sfm_data;
  if ( !Load( sfm_data, sSfM_Data_Filename, ESfM_Data( ALL ) ) )
  {
    OPENMVG_LOG_ERROR
              << "\nThe input SfM_Data file \"" << sSfM_Data_Filename << "\" cannot be read.";
    return EXIT_FAILURE;
  }

  if (floating_point_precision_digit < 6)
  {
    OPENMVG_LOG_WARNING << "\nThe precision will be set to 6\n";
    floating_point_precision_digit = 6;
  }

  if ( ! exportToColmap( sfm_data , sOutDir , floating_point_precision_digit) )
  {
    OPENMVG_LOG_ERROR << "There was an error during export of the file";
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
