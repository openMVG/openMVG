// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Romuald Perrot.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/cameras/Camera_Pinhole.hpp"
#include "openMVG/geometry/pose3.hpp"
#include "openMVG/geometry/frustum.hpp"
#include "openMVG/image/image_io.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_data_colorization.hpp"
#include "openMVG/stl/stl.hpp"

#include "config.h"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/progress/progress_display.hpp"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <iostream>
#include <fstream>

using namespace openMVG;
using namespace openMVG::sfm;
using namespace openMVG::geometry;
using namespace openMVG::image;
using namespace openMVG::cameras;

/**
* @brief Compute normal using the mean of the direction for each points
* @param sfm_data Input scene data
* @param vec_nor A vector of normal for each point (in same order as ColorizeTrack output)
*/
void ComputeNormals( const SfM_Data & sfm_data , std::vector<Vec3> & vec_nor )
{
  vec_nor.resize( sfm_data.GetLandmarks().size() );

  IndexT cpt = 0;
  for (Landmarks::const_iterator it = sfm_data.GetLandmarks().begin();
      it != sfm_data.GetLandmarks().end(); ++it, ++cpt)
    {
      // Set it's normal to 0,0,0
      vec_nor[ cpt ] = Vec3(0,0,0);

      const Landmark & landmark = it->second;
      const Observations & obs = landmark.obs;
      // Sum all rays from all observations
      for (Observations::const_iterator itOb = obs.begin(); itOb != obs.end(); ++itOb )
        {
          const IndexT viewId = itOb->first;

          // Compute ray of this obs
          const Observation & ob = itOb->second;
          const Vec2 & pos = ob.x;

          const View * view = sfm_data.GetViews().at(viewId).get();
          if (!sfm_data.IsPoseAndIntrinsicDefined(view))
          {
            continue;
          }
          // get it's intrinsic (if available)
          Intrinsics::const_iterator iterIntrinsic = sfm_data.GetIntrinsics().find(view->id_intrinsic);
          if (!isPinhole(iterIntrinsic->second->getType()))
          {
            continue;
          }

          const Pose3 pose = sfm_data.GetPoseOrDie(view);

          // Force a pinhole ? (TODO : need to work with another type of camera)
          const Pinhole_Intrinsic * cam = dynamic_cast<const Pinhole_Intrinsic*>(iterIntrinsic->second.get());

          // Compute vector from camera center to the point
          const Vec3 C = pose.center();
          const Mat3 Kinv = cam->K().inverse();
          const Mat3 Rt = pose.rotation().transpose();

          // TODO : good for educational but really inefficient ! (+C-C)
          const Vec3 dst = Rt * ( ( Kinv * Vec3( pos.x() , pos.y() , 1.0 ) ) ) + C;
          const Vec3 n = ( dst - C ).normalized();
          vec_nor[ cpt ] += n;
        }
        // Mean and normalize
        vec_nor[ cpt ] = vec_nor[ cpt ].normalized();
    }
}

/**
* @brief Export scene data to a file
* @param sfm_data Input scene data
* @param sOutFile Output file name
*/
bool ExportSceneData( const SfM_Data & sfm_data, const std::string & sOutFile )
{
  std::ofstream file( sOutFile );
  if ( ! file )
  {
    std::cerr << std::endl;
    std::cerr << "Could not create : " << sOutFile << std::endl;
    return false;
  }

  std::vector<Vec3> vec_3dPoints, vec_tracksColor , vec_normals;
  ColorizeTracks(sfm_data, vec_3dPoints, vec_tracksColor);
  ComputeNormals(sfm_data,vec_normals);

  // Model points
  file << "modelPos=[" << std::endl;
  for (size_t id_point = 0; id_point < vec_3dPoints.size(); ++id_point)
  {
    file << vec_3dPoints[id_point](0) << " , " << vec_3dPoints[id_point](1) << " , " << vec_3dPoints[id_point](2);
    if (id_point +1 != vec_3dPoints.size() )
    {
      file << " , ";
    }
  }
  file << "];" << std::endl;

  // Model normal
  file << "modelNor=[" << std::endl;
  for (size_t id_nor = 0; id_nor < vec_normals.size(); ++id_nor)
  {
    file << vec_normals[ id_nor ](0) << " , " << vec_normals[ id_nor ](1) << " , " << vec_normals[ id_nor ](2);
    if (id_nor + 1 != vec_normals.size() )
    {
      file << " , ";
    }
  }
  file << "];" << std::endl;

  // Model color
  file << "modelCol=[" << std::endl;
  for (size_t id_point = 0; id_point < vec_tracksColor.size();  ++id_point)
  {
    file << vec_tracksColor[id_point](0) << " , " << vec_tracksColor[id_point](1) << " , " << vec_tracksColor[id_point](2);
    if (id_point + 1 != vec_tracksColor.size() )
    {
      file << " , ";
    }
  }
  file << "];" << std::endl;

  // Camera positions
  file << "cameraPos=[" << std::endl;
  for (auto it = sfm_data.GetViews().begin(); it != sfm_data.GetViews().end();)
  {
    const View * view = sfm_data.GetViews().at(it->first).get();
    if (!sfm_data.IsPoseAndIntrinsicDefined(view))
    {
      ++it;
      continue;
    }
    // get it's intrinsic (if available)
    Intrinsics::const_iterator iterIntrinsic = sfm_data.GetIntrinsics().find(view->id_intrinsic);
    if (!isPinhole(iterIntrinsic->second->getType()))
    {
      ++it;
      continue;
    }

    const Pose3 pose = sfm_data.GetPoseOrDie(view);

    // Force a pinhole ? (TODO : need to work with another type of camera)
    const Pinhole_Intrinsic * cam = dynamic_cast<const Pinhole_Intrinsic*>(iterIntrinsic->second.get());
    if (cam == nullptr)
    {
      ++it;
      continue;
    }

    // Build frustum
    const Frustum f( cam->w() , cam->h() , cam->K() , pose.rotation() , pose.center() , 0.5 );

    const Vec3 pos = f.cones[0];
    file << "[" << pos[0] << "," << pos[1] << "," << pos[2] << "]";

    if ( ++it != sfm_data.GetViews().end() )
    {
      file << ",";
    }
  }

  file << "];" << std::endl;

  // Camera image planes
  file << "cameraImagePlanes=[" << std::endl;
  for (auto it = sfm_data.GetViews().begin(); it != sfm_data.GetViews().end();)
  {
    const View * view = sfm_data.GetViews().at(it->first).get();
    if (!sfm_data.IsPoseAndIntrinsicDefined(view))
    {
      ++it;
      continue;
    }
    // get it's intrinsic (if available)
    Intrinsics::const_iterator iterIntrinsic = sfm_data.GetIntrinsics().find(view->id_intrinsic);
    if (!isPinhole(iterIntrinsic->second->getType()))
    {
      ++it;
      continue;
    }

    const Pose3 pose = sfm_data.GetPoseOrDie(view);

    // Force a pinhole ? (TODO : need to work with another type of camera)
    const Pinhole_Intrinsic * cam = dynamic_cast<const Pinhole_Intrinsic*>(iterIntrinsic->second.get());
    if (cam == nullptr)
    {
      ++it;
      continue;
    }

    // Build frustum
    const Frustum f( cam->w() , cam->h() , cam->K() , pose.rotation() , pose.center() , 0.33 );

    const Vec3 & p1 = f.cones[1];
    const Vec3 & p2 = f.cones[2];
    const Vec3 & p3 = f.cones[3];
    const Vec3 & p4 = f.cones[4];

    file << "[";
    file << p1[0] << "," << p1[1] << "," << p1[2] << ",";
    file << p2[0] << "," << p2[1] << "," << p2[2] << ",";
    file << p3[0] << "," << p3[1] << "," << p3[2] << ",";
    file << p4[0] << "," << p4[1] << "," << p4[2] << "]";

    if ( ++it != sfm_data.GetViews().end() )
    {
      file << ",";
    }
  }
  file << "];" << std::endl;


  file.close();

  return true;
}

bool CopyFile( const std::string & inputFolder , const std::string & inputName ,
               const std::string & outputFolder ,
               const std::string & outputSubFolder , const std::string & outputName )
{
  std::string outputFilename;
  if (outputSubFolder.length() > 0 )
  {
    outputFilename = stlplus::create_filespec( stlplus::folder_append_separator( outputFolder ) + outputSubFolder , outputName );
  }
  else
  {
    outputFilename = stlplus::create_filespec( outputFolder , outputName );
  }

  const std::string inputFilename1 = stlplus::create_filespec( stlplus::folder_append_separator( OPENMVG_SFM_WEBGL_SRC_PATH ) + inputFolder , inputName );
  const std::string inputFilename2 = stlplus::create_filespec( stlplus::folder_append_separator( OPENMVG_SFM_WEBGL_INSTALL_SHARED_PATH ) +  inputFolder , inputName );

  // Prefer Build path
  std::string usablePath = inputFilename1;
  if ( ! stlplus::file_exists( inputFilename1 ) )
  {
    usablePath = inputFilename2;
  }

  if ( ! stlplus::file_copy( usablePath , outputFilename ) )
  {
    std::cerr << "Error copying " << inputName << " file" << std::endl;
    return false;
  }
  return true;
}

/**
* Stupid template trick to get the size of a static c array
*/
template < typename T , int N>
int StaticArraySize( T (&)[N] )
{
  return N;
}

/**
* @brief Export a sfm scene to a webgl app
* @param sSfM_Data_Filename Filename of the sfm scene
* @param sOutDir Output directory
* @note If output directory no present, create it
* @retval true if everything was correct
* @retval false if there was an error
*/
bool ExportToWebGL( const std::string & sSfM_Data_Filename , const std::string & sOutDir )
{
  // Load input scene
  SfM_Data sfm_data;
  if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(ALL)))
  {
    std::cerr << std::endl
      << "The input SfM_Data file \""<< sSfM_Data_Filename << "\" cannot be read." << std::endl;
    return false;
  }


  // Check output directory
  // OUT_DIR /
  //         / common
  //         / style
  //         / model
  bool bOk = true;
  if (!stlplus::is_folder(sOutDir))
  {
    stlplus::folder_create(sOutDir);
    bOk = stlplus::is_folder(sOutDir);
  }
  stlplus::folder_create( stlplus::folder_append_separator(sOutDir) + "common" );
  stlplus::folder_create( stlplus::folder_append_separator(sOutDir) + "style" );
  stlplus::folder_create( stlplus::folder_append_separator(sOutDir) + "model" );

  if (! stlplus::is_folder( stlplus::folder_append_separator(sOutDir) + "common") ||
      ! stlplus::is_folder( stlplus::folder_append_separator(sOutDir) + "style") ||
      ! stlplus::is_folder( stlplus::folder_append_separator(sOutDir) + "model" ) )
  {
    bOk = false;
  }
  if ( ! bOk )
  {
    std::cerr << "Could not create directory structure : \"" << sOutDir << "\"" << std::endl;
    return false;
  }


  // 1 - Export scene specific data
  const std::string modelFilename = stlplus::create_filespec(
    stlplus::folder_append_separator( sOutDir ) + "model" , "model.js" );
    std::cerr << "Modelfilename:" << modelFilename << std::endl;
  ExportSceneData( sfm_data , modelFilename );

  const std::string files_to_copy[] =
  {
    "common" , "index.html" , sOutDir , "" , "index.html" , // 1
    "common" , "style.css" , sOutDir , "style" , "style.css" , // 2
    "common" , "main.js" , sOutDir , "common" , "main.js" , // 3
    "common" , "quaternion.js" , sOutDir , "common" , "quaternion.js" , // 4
    "common" , "dual_quaternion.js" , sOutDir , "common" , "dual_quaternion.js" , // 5
    "common" , "common.js" , sOutDir , "common" , "common.js" , // 6
    "common" , "matrix4.js" , sOutDir , "common" , "matrix4.js" , // 7
    "common" , "matrix3.js" , sOutDir , "common" , "matrix3.js" , // 8
    "common" , "plane.js" , sOutDir , "common" , "plane.js" , // 9
    "common" , "camera.js" , sOutDir , "common" , "camera.js" , // 10
    "common" , "shader.js" , sOutDir , "common" , "shader.js" , // 11
    "common" , "trackball.js" , sOutDir , "common" , "trackball.js" , // 12
    "common" , "render_context.js" , sOutDir , "common" , "render_context.js" , // 13
    "common" , "vector.js" , sOutDir , "common" , "vector.js" , // 14
    "common" , "point_cloud.js" , sOutDir , "common" , "point_cloud.js" , // 15
    "common" , "camera_gizmo.js" , sOutDir , "common" , "camera_gizmo.js" // 16
  };

  // Copy the js files
  for (int id_file = 0; id_file < StaticArraySize( files_to_copy ) / 5; ++id_file)
  {
    const std::string & input_folder = files_to_copy[ 5 * id_file ];
    const std::string & input_name   = files_to_copy[ 5 * id_file + 1 ];
    const std::string & output_folder = files_to_copy[ 5 * id_file + 2 ];
    const std::string & output_subfolder = files_to_copy[ 5 * id_file + 3 ];
    const std::string & output_name = files_to_copy[ 5 * id_file + 4 ];

    if ( ! CopyFile( input_folder , input_name , output_folder , output_subfolder , output_name ) )
    {
      return false;
    }
  }

  return true;
}

int main( int argc , char ** argv )
{
  // -i sfm_data.bin -o directory
  CmdLine cmd;
  std::string sSfM_Data_Filename;
  std::string sOutDir = "";

  cmd.add( make_option('i', sSfM_Data_Filename, "sfmdata") );
  cmd.add( make_option('o', sOutDir, "outdir") );


  try
  {
      if (argc == 1)
      {
        throw std::string("Invalid command line parameter.");
      }
      cmd.process(argc, argv);
  }
  catch (const std::string& s)
  {
      std::cerr << "Usage: " << argv[0] << '\n'
      << "[-i|--sfmdata] filename, the SfM_Data file to convert\n"
      << "[-o|--outdir path]\n"
      << std::endl;

      std::cerr << s << std::endl;
      exit( EXIT_FAILURE );
  }

  std::cout << "Using : " << std::endl;
  std::cout << "Input file : " << sSfM_Data_Filename << std::endl;
  std::cout << "Output dir : " << sOutDir << std::endl;

  if ( ! ExportToWebGL( sSfM_Data_Filename , sOutDir ) )
  {
    std::cerr << "There was an error during export" << std::endl;
    exit( EXIT_FAILURE );
  }


  return EXIT_SUCCESS;
}
