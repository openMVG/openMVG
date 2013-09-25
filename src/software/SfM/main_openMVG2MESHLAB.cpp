#include "software/SfMViewer/document.h"

#include "openMVG/image/image.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "openMVG/numeric/numeric.h"

#include <fstream> 


int main(int argc, char **argv)
{
  CmdLine cmd;

  std::string sSfM_OutputPath;
  std::string sPly = "";
  std::string sOutDir = "";

  cmd.add( make_option('s', sSfM_OutputPath, "SfMPath") );
  cmd.add( make_option('p', sPly, "ply") );
  cmd.add( make_option('o', sOutDir, "outdir") );

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch(const std::string& s) {
      std::cerr << "Usage: " << argv[0] << ' '
      << "[-s|--SfMPath] "
      << "[-p|--ply path] "
      << "[-o|--outdir path] "
      << std::endl;

      std::cerr << s << std::endl;
      return EXIT_FAILURE;
  }

  std::cout << " You called : " <<std::endl
            << argv[0] << std::endl
            << "--SfMPath " << sSfM_OutputPath << std::endl
            << "--ply " << sPly << std::endl
            << "--outdir " << sOutDir << std::endl;

  if (!stlplus::folder_exists(sSfM_OutputPath) )  {
    std::cerr << "\nSfM directory doesn't exist" << std::endl;
    return EXIT_FAILURE;
  } 
  
  if (!stlplus::file_exists(sPly) )  {
    std::cerr << "\nPly file doesn't exist" << std::endl;
    return EXIT_FAILURE;
  } 

  // Create output dir
  if (!stlplus::folder_exists(sOutDir))
    stlplus::folder_create( sOutDir );

  Document m_doc;
  std::cout << "\n Open the directory : \n" << sSfM_OutputPath << std::endl;

  //Read SfM output directory
  if (!m_doc.load(sSfM_OutputPath))
  {
    std::cout << "Impossible to open the openMVG SfM project." << std::endl;
    return EXIT_FAILURE;
  }
  
  std::ofstream outfile( stlplus::create_filespec( sOutDir, "sceneMeshlab", "mlp" ).c_str() );
  
  //Init mlp file  
  outfile << "<!DOCTYPE MeshLabDocument>" << std::endl 
    << "<MeshLabProject>" << std::endl 
    << " <MeshGroup>" << std::endl 
    << "  <MLMesh label=\"" << sPly << "\" filename=\"" << sPly << "\">" << std::endl 
    << "   <MLMatrix44>" << std::endl 
    << "1 0 0 0 " << std::endl 
    << "0 1 0 0 " << std::endl 
    << "0 0 1 0 " << std::endl 
    << "0 0 0 1 " << std::endl 
    << "</MLMatrix44>" << std::endl 
    << "  </MLMesh>" << std::endl 
    << " </MeshGroup>" << std::endl; 
  
  outfile <<  " <RasterGroup>" << std::endl;
  
  std::map<size_t, PinholeCamera >::const_iterator iterCamera = m_doc._map_camera.begin();
  std::map<size_t, std::pair<size_t,size_t> >::const_iterator iterSize = m_doc._map_imageSize.begin();
  std::vector<std::string>::const_iterator iterName = m_doc._vec_imageNames.begin();
  for ( ;
        iterCamera != m_doc._map_camera.end();
        iterCamera++,
        iterSize++,
        iterName++ )
  {
    PinholeCamera camera = iterCamera->second;
    Mat34 P = camera._P;
    for ( int i = 1; i < 3 ; i++)
      for ( int j = 0; j < 4; j++)
        P(i, j) *= -1;

    Mat3 R, K;
    Vec3 t;
    KRt_From_P( P, &K, &R, &t);

    Vec3 optical_center = R.transpose() * t;

    outfile << "  <MLRaster label=\"" << *iterName << "\">" << std::endl
      << "   <VCGCamera TranslationVector=\""
      << optical_center[0] << " " 
      << optical_center[1] << " " 
      << optical_center[2] << " " 
      << " 1 \"" 
      << " LensDistortion=\"0 0\""
      << " ViewportPx=\"" << iterSize->second.first << " " << iterSize->second.second << "\"" 
      << " PixelSizeMm=\"" << 1  << " " << 1 << "\""  
      << " CenterPx=\"" << iterSize->second.first / 2.0 << " " << iterSize->second.second / 2.0 << "\""
      << " FocalMm=\"" << (double)K(0, 0 )  << "\"" 
      << " RotationMatrix=\""
      << R(0, 0) << " " << R(0, 1) << " " << R(0, 2) << " 0 "
      << R(1, 0) << " " << R(1, 1) << " " << R(1, 2) << " 0 "
      << R(2, 0) << " " << R(2, 1) << " " << R(2, 2) << " 0 " 
      << "0 0 0 1 \"/>"  << std::endl;
    std::string soffsetImagePath =  stlplus::create_filespec( sSfM_OutputPath, "imagesOffset" );
    if ( stlplus::folder_exists( soffsetImagePath ) )
      outfile << "   <Plane semantic=\"\" fileName=\"" << stlplus::create_filespec( soffsetImagePath, 
        stlplus::basename_part( *iterName ) + "_of", 
        stlplus::extension_part( *iterName ) ) << "\"/> "<< std::endl;
    else
      outfile << "   <Plane semantic=\"\" fileName=\"" << stlplus::create_filespec( stlplus::create_filespec( sSfM_OutputPath, "images" ), 
        *iterName ) << "\"/> "<< std::endl;      
      outfile << "  </MLRaster>" << std::endl;
  }
  
  outfile << "   </RasterGroup>" << std::endl
    << "</MeshLabProject>" << std::endl;

  outfile.close();

  return EXIT_SUCCESS;
}