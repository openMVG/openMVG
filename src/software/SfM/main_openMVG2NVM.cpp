#include "openMVG/sfm/sfm.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/progress/progress.hpp"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

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
                    const std::string sOutDirectory ,
                    const std::string filename )
{
  const string sOutViewsDirectory2 = "views" ;
  const string sOutViewsDirectory = stlplus::folder_append_separator( sOutDirectory ) + "views";
  if ( !stlplus::folder_exists( sOutViewsDirectory ) )
  {
    cout << "\033[1;31mCreating directory:  " << sOutViewsDirectory << "\033[0m\n";
    stlplus::folder_create( sOutViewsDirectory );
  }

  // Header
  std::ofstream file( filename ) ;

  if( ! file )
  {
    std::cerr << "Cannot write file" << filename << std::endl ;
    return false ;
  }
  file << "NVM_V3" << std::endl ;

  // we reindex the poses to ensure a contiguous pose list.
  Hash_Map<IndexT, IndexT> map_viewIdToContiguous;
  int nb_cam = 0 ;
  for( Views::const_iterator iter = sfm_data.GetViews().begin();
       iter != sfm_data.GetViews().end(); ++iter )
  {
    const View * view = iter->second.get();
    if ( !sfm_data.IsPoseAndIntrinsicDefined( view ) )
    {
      continue;
    }

    nb_cam ++ ;
    map_viewIdToContiguous.insert( std::make_pair( view->id_view, map_viewIdToContiguous.size() ) );
  }

  // Number of cameras
  // For each camera : File_name Focal Qw Qx Qy Qz Cx Cy Cz D0 D1 0

  file << nb_cam << std::endl ;
  std::string sOutViewIteratorDirectory;
  std::string sOutViewIteratorDirectory2 ;
  Image<RGBColor> image, image_ud ;
  C_Progress_display my_progress_bar( sfm_data.GetViews().size() );
  for( Views::const_iterator iter = sfm_data.GetViews().begin();
       iter != sfm_data.GetViews().end(); ++iter , ++my_progress_bar )
  {
    const View * view = iter->second.get();

    if ( !sfm_data.IsPoseAndIntrinsicDefined( view ) )
    {
      continue;
    }

    std::ostringstream padding;
    padding << std::setw( 4 ) << std::setfill( '0' ) << view->id_view;
    sOutViewIteratorDirectory = stlplus::folder_append_separator( sOutViewsDirectory ) + "view_" + padding.str() ;
    sOutViewIteratorDirectory2 = stlplus::folder_append_separator( sOutViewsDirectory2 ) + "view_" + padding.str() ;

    const std::string srcImage = stlplus::create_filespec( sfm_data.s_root_path, view->s_Img_path );
    const std::string dstImage =
      stlplus::create_filespec( stlplus::folder_append_separator( sOutViewIteratorDirectory ), "undistorted", "png" );
    const std::string dstImage2 =
      stlplus::create_filespec( stlplus::folder_append_separator( sOutViewIteratorDirectory2 ), "undistorted", "png" );

    // Create output dir if not present
    if ( !stlplus::folder_exists( sOutViewIteratorDirectory ) )
    {
      stlplus::folder_create( sOutViewIteratorDirectory );
    }



    Intrinsics::const_iterator iterIntrinsic = sfm_data.GetIntrinsics().find( view->id_intrinsic );
    const IntrinsicBase * cam = iterIntrinsic->second.get();
    const Pinhole_Intrinsic * pinhole_cam = static_cast<const Pinhole_Intrinsic *>( cam );
    const double flen = pinhole_cam->focal() ;
    const Pose3 pose = sfm_data.GetPoseOrDie( view );
    const Mat3 rotation = pose.rotation();
    const Vec3 center = pose.center();

    const double Cx = center[0] ;
    const double Cy = center[1] ;
    const double Cz = center[2] ;
    Eigen::Quaterniond q( rotation ) ;
    const double Qx = q.x() ;
    const double Qy = q.y() ;
    const double Qz = q.z() ;
    const double Qw = q.w() ;
    double d0 = 0.0 ;
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


    file << dstImage2 << " " <<
         flen << " " <<

         Qw << " " <<
         Qx << " " <<
         Qy << " " <<
         Qz << " " <<

         Cx << " " <<
         Cy << " " <<
         Cz << " " <<

         d0 << " " <<
         0 << std::endl ;
  }

  // Now exports points
  // Number of points
  // For each points : X Y Z R G B Nm [ mesurements ]
  // mesurements : Img_idx Feat_idx X Y
  const Landmarks & landmarks = sfm_data.GetLandmarks();
  const size_t featureCount = std::distance( landmarks.begin(), landmarks.end() );
  file << featureCount << std::endl ;
  for ( Landmarks::const_iterator iterLandmarks = landmarks.begin(); iterLandmarks != landmarks.end(); ++iterLandmarks )
  {
    const Vec3 exportPoint = iterLandmarks->second.X;
    file << exportPoint.x() << " " << exportPoint.y() << " " << exportPoint.z() << " ";
    file << 250 << " " << 100 << " " << 150 << " ";  // Write arbitrary RGB color, see above note

    // Tally set of feature observations
    const Observations & obs = iterLandmarks->second.obs;
    const size_t featureCount = std::distance( obs.begin(), obs.end() );
    file << featureCount;

    for ( Observations::const_iterator itObs = obs.begin(); itObs != obs.end(); ++itObs )
    {
      const IndexT viewId = map_viewIdToContiguous[itObs->first] ;
      const IndexT featId = itObs->second.id_feat;
      const Observation & ob = itObs->second;

      file << " " << viewId << " " << featId << " " << ob.x( 0 ) << " " << ob.x( 1 ) << " " ;
    }
    file << "\n";
  }



  return true ;
}

/**
* @brief Main function used to export a NVM file
* @param sfm_data Structure from Motion file to export
* @param sOutDirectory Output directory
*/
bool exportToNVM( const SfM_Data & sfm_data , const std::string & sOutDirectory  )
{
  // Create output directory
  bool bOk = false ;
  if ( !stlplus::is_folder( sOutDirectory ) )
  {
    cout << "\033[1;31mCreating directory:  " << sOutDirectory << "\033[0m\n";
    stlplus::folder_create( sOutDirectory );
    bOk = stlplus::is_folder( sOutDirectory );
  }
  else
  {
    bOk = true ;
  }

  if ( !bOk )
  {
    std::cerr << "Cannot access one of the desired output directories" << std::endl;
    return false;
  }
  const std::string sFilename = stlplus::create_filespec( sOutDirectory , "scene.nvm" ) ;
  if( ! CreateNVMFile( sfm_data , sOutDirectory , sFilename ) )
  {
    std::cerr << "There was an error exporting project" << std::endl ;
    return false ;
  }
  return true ;

}

int main( int argc , char ** argv )
{
  CmdLine cmd;
  std::string sSfM_Data_Filename;
  std::string sOutDir = "";
  cmd.add( make_option( 'i', sSfM_Data_Filename, "sfmdata" ) );
  cmd.add( make_option( 'o', sOutDir, "outdir" ) );
  std::cout << "Note:  this program writes output in NVM file format.\n";

  try
  {
    if ( argc == 1 )
    {
      throw std::string( "Invalid command line parameter." );
    }
    cmd.process( argc, argv );
  }
  catch( const std::string& s )
  {
    std::cerr << "Usage: " << argv[0] << '\n'
              << "[-i|--sfmdata] filename, the SfM_Data file to convert\n"
              << "[-o|--outdir] path\n"
              << std::endl;

    std::cerr << s << std::endl;
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
    std::cerr << std::endl
              << "The input SfM_Data file \"" << sSfM_Data_Filename << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }

  if( ! exportToNVM( sfm_data , sOutDir ) )
  {
    std::cerr << "There was an error during export of the file" << std::endl;
    exit( EXIT_FAILURE ) ;
  }


  return EXIT_SUCCESS ;
}