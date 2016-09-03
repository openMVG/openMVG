#include <iostream>

#include "third_party/cmdLine/cmdLine.h"
#include "openMVG/sfm/sfm.hpp"
#include "openMVG/stl/stl.hpp"


#include "config.h"

using namespace openMVG;
using namespace openMVG::sfm;
using namespace openMVG::image;


/// Find the color of the SfM_Data Landmarks/structure
bool ColorizeTracks(
  const SfM_Data & sfm_data,
  std::vector<Vec3> & vec_3dPoints,
  std::vector<Vec3> & vec_tracksColor)
{
  // Colorize each track
  //  Start with the most representative image
  //    and iterate to provide a color to each 3D point

  {
    C_Progress_display my_progress_bar(sfm_data.GetLandmarks().size(),
                                       std::cout,
                                       "\nCompute scene structure color\n");

    vec_tracksColor.resize(sfm_data.GetLandmarks().size());
    vec_3dPoints.resize(sfm_data.GetLandmarks().size());

    //Build a list of contiguous index for the trackIds
    std::map<IndexT, IndexT> trackIds_to_contiguousIndexes;
    IndexT cpt = 0;
    for (Landmarks::const_iterator it = sfm_data.GetLandmarks().begin();
      it != sfm_data.GetLandmarks().end(); ++it, ++cpt)
    {
      trackIds_to_contiguousIndexes[it->first] = cpt;
      vec_3dPoints[cpt] = it->second.X;
    }

    // The track list that will be colored (point removed during the process)
    std::set<IndexT> remainingTrackToColor;
    std::transform(sfm_data.GetLandmarks().begin(), sfm_data.GetLandmarks().end(),
      std::inserter(remainingTrackToColor, remainingTrackToColor.begin()),
      stl::RetrieveKey());

    while( !remainingTrackToColor.empty() )
    {
      // Find the most representative image (for the remaining 3D points)
      //  a. Count the number of observation per view for each 3Dpoint Index
      //  b. Sort to find the most representative view index

      std::map<IndexT, IndexT> map_IndexCardinal; // ViewId, Cardinal
      for (std::set<IndexT>::const_iterator
        iterT = remainingTrackToColor.begin();
        iterT != remainingTrackToColor.end();
        ++iterT)
      {
        const size_t trackId = *iterT;
        const Observations & obs = sfm_data.GetLandmarks().at(trackId).obs;
        for( Observations::const_iterator iterObs = obs.begin();
          iterObs != obs.end(); ++iterObs)
        {
          const size_t viewId = iterObs->first;
          if (map_IndexCardinal.find(viewId) == map_IndexCardinal.end())
            map_IndexCardinal[viewId] = 1;
          else
            ++map_IndexCardinal[viewId];
        }
      }

      // Find the View index that is the most represented
      std::vector<IndexT> vec_cardinal;
      std::transform(map_IndexCardinal.begin(),
        map_IndexCardinal.end(),
        std::back_inserter(vec_cardinal),
        stl::RetrieveValue());
      using namespace stl::indexed_sort;
      std::vector< sort_index_packet_descend< IndexT, IndexT> > packet_vec(vec_cardinal.size());
      sort_index_helper(packet_vec, &vec_cardinal[0], 1);

      // First image index with the most of occurence
      std::map<IndexT, IndexT>::const_iterator iterTT = map_IndexCardinal.begin();
      std::advance(iterTT, packet_vec[0].index);
      const size_t view_index = iterTT->first;
      const View * view = sfm_data.GetViews().at(view_index).get();
      const std::string sView_filename = stlplus::create_filespec(sfm_data.s_root_path,
        view->s_Img_path);
      Image<RGBColor> image_rgb;
      Image<unsigned char> image_gray;
      const bool b_rgb_image = ReadImage(sView_filename.c_str(), &image_rgb);
      if (!b_rgb_image) //try Gray level
      {
        const bool b_gray_image = ReadImage(sView_filename.c_str(), &image_gray);
        if (!b_gray_image)
        {
          std::cerr << "Cannot open provided the image." << std::endl;
          return false;
        }
      }

      // Iterate through the remaining track to color
      // - look if the current view is present to color the track
      std::set<IndexT> set_toRemove;
      for (std::set<IndexT>::const_iterator
        iterT = remainingTrackToColor.begin();
        iterT != remainingTrackToColor.end();
        ++iterT)
      {
        const size_t trackId = *iterT;
        const Observations & obs = sfm_data.GetLandmarks().at(trackId).obs;
        Observations::const_iterator it = obs.find(view_index);

        if (it != obs.end())
        {
          // Color the track
          const Vec2 & pt = it->second.x;
          const RGBColor color = b_rgb_image ? image_rgb(pt.y(), pt.x()) : RGBColor(image_gray(pt.y(), pt.x()));

          vec_tracksColor[ trackIds_to_contiguousIndexes[trackId] ] = Vec3(color.r(), color.g(), color.b());
          set_toRemove.insert(trackId);
          ++my_progress_bar;
        }
      }
      // Remove colored track
      for (std::set<IndexT>::const_iterator iter = set_toRemove.begin();
        iter != set_toRemove.end(); ++iter)
      {
        remainingTrackToColor.erase(*iter);
      }
    }
  }
  return true;
}

/**
* @brief Export scene data to a file 
* @param sfm_data Input scene data 
* @param sOutFile Output file name 
*/
bool ExportSceneData( const SfM_Data & sfm_data, const std::string & sOutFile )
{
  std::ofstream file( sOutFile ) ;
  if( ! file )
  {
    std::cerr << std::endl ;
    std::cerr << "Could not create : " << sOutFile << std::endl ; 
    return false ; 
  }

  std::vector<Vec3> vec_3dPoints, vec_tracksColor ;
  ColorizeTracks(sfm_data, vec_3dPoints, vec_tracksColor) ;

  // 3d points 
  file << "modelPos=[" << std::endl ; 
  for( size_t id_point = 0 ; id_point < vec_3dPoints.size() ; ++id_point )
  {
    file << vec_3dPoints[id_point](0) << " , " << vec_3dPoints[id_point](1) << " , " << vec_3dPoints[id_point](2) ;
    if( id_point +1 != vec_3dPoints.size() )
    {
      file << " , " ;
    }
  }
  file << "];" << std::endl ; 

  file << "modelCol=[" << std::endl ; 
  for( size_t id_point = 0 ; id_point < vec_tracksColor.size() ;  ++id_point )
  {
    file << vec_tracksColor[id_point](0) << " , " << vec_tracksColor[id_point](1) << " , " << vec_tracksColor[id_point](2) ;
    if( id_point + 1 != vec_tracksColor.size() )
    {
      file << " , " ;
    }
  }
  file << "];" << std::endl ; 

  file.close() ; 
  
  return true ; 
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
  stlplus::folder_create( stlplus::folder_append_separator(sOutDir) + "common" ) ;
  stlplus::folder_create( stlplus::folder_append_separator(sOutDir) + "style" ) ;
  stlplus::folder_create( stlplus::folder_append_separator(sOutDir) + "model" ) ; 
  
  if( ! stlplus::is_folder( stlplus::folder_append_separator(sOutDir) + "common") ||
      ! stlplus::is_folder( stlplus::folder_append_separator(sOutDir) + "style") ||
      ! stlplus::is_folder( stlplus::folder_append_separator(sOutDir) + "model" ) )
  {
    bOk = false ; 
  }
  if( ! bOk )
  {
    std::cerr << "Could not create directory structure : \"" << sOutDir << "\"" << std::endl ;
    return false ; 
  }


  // 1 - Export scene specific data  
  const std::string modelFilename = stlplus::create_filespec( 
    stlplus::folder_append_separator( sOutDir ) + "model" , "model.js" ) ;
    std::cerr << "Modelfilename:" << modelFilename << std::endl ; 
  ExportSceneData( sfm_data , modelFilename ) ; 


  // 2 - Copy common files 
  // - index 
  const std::string outputIndexFilename = stlplus::create_filespec( sOutDir , "index.html" ) ; 
  const std::string inputIndexFilename = stlplus::create_filespec( 
    stlplus::folder_append_separator( OPENMVG_SFM_WEBGL_PATH ) + "common" , "index.html" ) ;
  if( ! stlplus::file_copy( inputIndexFilename , outputIndexFilename ) )
  {
    std::cerr << "Error copying index.html file" << std::endl ; 
    return false ; 
  }

  // - style 
  const std::string outputStyleFilename = stlplus::create_filespec( 
    stlplus::folder_append_separator( sOutDir ) + "style" , "style.css" ) ;
  const std::string inputStyleFilename = stlplus::create_filespec( 
    stlplus::folder_append_separator( OPENMVG_SFM_WEBGL_PATH ) + "common" , "style.css" ) ;
  if( ! stlplus::file_copy( inputStyleFilename , outputStyleFilename ) )
  {
    std::cerr << "Error copying style/style.css file" << std::endl ; 
    return false ; 
  }

  // - main  
  const std::string outputMainFilename = stlplus::create_filespec( 
    stlplus::folder_append_separator( sOutDir ) + "common" , "main.js" ) ;
  const std::string inputMainFilaname = stlplus::create_filespec( 
    stlplus::folder_append_separator( OPENMVG_SFM_WEBGL_PATH ) + "common" , "main.js" ) ;
  if( ! stlplus::file_copy( inputMainFilaname , outputMainFilename ) )
  {
    std::cerr << "Error copying common/main.js file" << std::endl ; 
    return false ; 
  }
  // - quaternion 
  const std::string outputQuaternionFilename = stlplus::create_filespec( 
    stlplus::folder_append_separator( sOutDir ) + "common" , "quaternion.js" ) ;
  const std::string inputQuaternionFilaname = stlplus::create_filespec( 
    stlplus::folder_append_separator( OPENMVG_SFM_WEBGL_PATH ) + "common" , "quaternion.js" ) ;
  if( ! stlplus::file_copy( inputQuaternionFilaname , outputQuaternionFilename ) )
  {
    std::cerr << "Error copying common/quaternion.js file" << std::endl ; 
    return false ; 
  }
  // - dual quaternion 
  const std::string outputDualQuaternionFilename = stlplus::create_filespec( 
    stlplus::folder_append_separator( sOutDir ) + "common" , "dual_quaternion.js" ) ;
  const std::string inputDualQuaternionFilaname = stlplus::create_filespec( 
    stlplus::folder_append_separator( OPENMVG_SFM_WEBGL_PATH ) + "common" , "dual_quaternion.js" ) ;
  if( ! stlplus::file_copy( inputDualQuaternionFilaname , outputDualQuaternionFilename ) )
  {
    std::cerr << "Error copying common/dual_quaternion.js file" << std::endl ; 
    return false ; 
  }
  // - common  
  const std::string outputCommonFilename = stlplus::create_filespec( 
    stlplus::folder_append_separator( sOutDir ) + "common" , "common.js" ) ;
  const std::string inputCommonFilaname = stlplus::create_filespec( 
    stlplus::folder_append_separator( OPENMVG_SFM_WEBGL_PATH ) + "common" , "common.js" ) ;
  if( ! stlplus::file_copy( inputCommonFilaname , outputCommonFilename ) )
  {
    std::cerr << "Error copying common/common.js file" << std::endl ; 
    return false ; 
  }
  // - matrix3 
  const std::string outputMatrix3Filename = stlplus::create_filespec( 
    stlplus::folder_append_separator( sOutDir ) + "common" , "matrix3.js" ) ;
  const std::string inputMatrix3Filaname = stlplus::create_filespec( 
    stlplus::folder_append_separator( OPENMVG_SFM_WEBGL_PATH ) + "common" , "matrix3.js" ) ;
  if( ! stlplus::file_copy( inputMatrix3Filaname , outputMatrix3Filename ) )
  {
    std::cerr << "Error copying common/matrix3.js file" << std::endl ; 
    return false ; 
  }
  // - matrix4 
  const std::string outputMatrix4Filename = stlplus::create_filespec( 
    stlplus::folder_append_separator( sOutDir ) + "common" , "matrix4.js" ) ;
  const std::string inputMatrix4Filaname = stlplus::create_filespec( 
    stlplus::folder_append_separator( OPENMVG_SFM_WEBGL_PATH ) + "common" , "matrix4.js" ) ;
  if( ! stlplus::file_copy( inputMatrix4Filaname , outputMatrix4Filename ) )
  {
    std::cerr << "Error copying common/matrix4.js file" << std::endl ; 
    return false ; 
  }
  // - plane 
  const std::string outputPlaneFilename = stlplus::create_filespec( 
  stlplus::folder_append_separator( sOutDir ) + "common" , "plane.js" ) ;
  const std::string inputPlaneFilaname = stlplus::create_filespec( 
    stlplus::folder_append_separator( OPENMVG_SFM_WEBGL_PATH ) + "common" , "plane.js" ) ;
  if( ! stlplus::file_copy( inputPlaneFilaname , outputPlaneFilename ) )
  {
    std::cerr << "Error copying common/plane.js file" << std::endl ; 
    return false ; 
  }
  // - camera 
  const std::string outputCameraFilename = stlplus::create_filespec( 
    stlplus::folder_append_separator( sOutDir ) + "common" , "camera.js" ) ;
  const std::string inputCameraFilaname = stlplus::create_filespec( 
    stlplus::folder_append_separator( OPENMVG_SFM_WEBGL_PATH ) + "common" , "camera.js" ) ;
  if( ! stlplus::file_copy( inputCameraFilaname , outputCameraFilename ) )
  {
    std::cerr << "Error copying common/camera.js file" << std::endl ; 
    return false ; 
  }
  // - shader 
  const std::string outputShaderFilename = stlplus::create_filespec( 
    stlplus::folder_append_separator( sOutDir ) + "common" , "shader.js" ) ;
  const std::string inputShaderFilaname = stlplus::create_filespec( 
    stlplus::folder_append_separator( OPENMVG_SFM_WEBGL_PATH ) + "common" , "shader.js" ) ;
  if( ! stlplus::file_copy( inputShaderFilaname , outputShaderFilename ) )
  {
    std::cerr << "Error copying common/shader.js file" << std::endl ; 
    return false ; 
  }
  // - trackball 
  const std::string outputTrackballFilename = stlplus::create_filespec( 
    stlplus::folder_append_separator( sOutDir ) + "common" , "trackball.js" ) ;
  const std::string inputTrackballFilaname = stlplus::create_filespec( 
    stlplus::folder_append_separator( OPENMVG_SFM_WEBGL_PATH ) + "common" , "trackball.js" ) ;
  if( ! stlplus::file_copy( inputTrackballFilaname , outputTrackballFilename ) )
  {
    std::cerr << "Error copying common/trackball.js file" << std::endl ; 
    return false ; 
  }
  // - render_context 
  const std::string outputRenderContextFilename = stlplus::create_filespec( 
    stlplus::folder_append_separator( sOutDir ) + "common" , "render_context.js" ) ;
  const std::string inputRenderContextFilaname = stlplus::create_filespec( 
    stlplus::folder_append_separator( OPENMVG_SFM_WEBGL_PATH ) + "common" , "render_context.js" ) ;
  if( ! stlplus::file_copy( inputRenderContextFilaname , outputRenderContextFilename ) )
  {
    std::cerr << "Error copying common/render_context.js file" << std::endl ; 
    return false ; 
  }
  // - vector 
  const std::string outputVectorFilename = stlplus::create_filespec( 
    stlplus::folder_append_separator( sOutDir ) + "common" , "vector.js" ) ;
  const std::string inputVectorFilaname = stlplus::create_filespec( 
    stlplus::folder_append_separator( OPENMVG_SFM_WEBGL_PATH ) + "common" , "vector.js" ) ;
  if( ! stlplus::file_copy( inputVectorFilaname , outputVectorFilename ) )
  {
    std::cerr << "Error copying common/vector.js file" << std::endl ; 
    return false ; 
  }
  // - point_cloud
  const std::string outputPointCloudFilename = stlplus::create_filespec( 
    stlplus::folder_append_separator( sOutDir ) + "common" , "point_cloud.js" ) ;
  const std::string inputPointCloudFilaname = stlplus::create_filespec( 
    stlplus::folder_append_separator( OPENMVG_SFM_WEBGL_PATH ) + "common" , "point_cloud.js" ) ;
  if( ! stlplus::file_copy( inputPointCloudFilaname , outputPointCloudFilename ) )
  {
    std::cerr << "Error copying common/point_cloud.js file" << std::endl ; 
    return false ; 
  }
  

  return true ; 
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
  catch(const std::string& s) 
  {
      std::cerr << "Usage: " << argv[0] << '\n'
      << "[-i|--sfmdata] filename, the SfM_Data file to convert\n"
      << "[-o|--outdir path]\n"
      << std::endl;

      std::cerr << s << std::endl;
      exit( EXIT_FAILURE ) ;
  }

  std::cout << "Using : " << std::endl ; 
  std::cout << "Input file : " << sSfM_Data_Filename << std::endl ; 
  std::cout << "Output dir : " << sOutDir << std::endl ; 

  if( ! ExportToWebGL( sSfM_Data_Filename , sOutDir ) )
  {
    std::cerr << "There was an error during export" << std::endl ; 
    exit( EXIT_FAILURE ) ;  
  }
  

  return EXIT_SUCCESS ;
}