// Copyright (c) 2014 Romuald Perrot, Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/image/image.hpp"
#include "openMVG/features/akaze/AKAZE.hpp"

#include "third_party/vectorGraphics/svgDrawer.hpp"
using namespace svg;


#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include "openMVG/system/timer.hpp"

#include <iostream>
#include <sstream>

using namespace openMVG;

void usage( const std::string & appName )
{
  std::cerr << "usage : " << std::endl ;
  std::cerr << "  " << appName << std::endl ;
  std::cerr << "[-i|--input imageFileName]" << std::endl ;
  std::cerr << "[-o|--output outputFileName]" << std::endl ;
  std::cerr << "[Optional]" << std::endl ;
  std::cerr << "[-p|--nb-octave 4]" << std::endl ;
  std::cerr << "[-q|--nb-slice 4]" << std::endl ;

  exit( EXIT_FAILURE ) ;
}

int main( int argc , char ** argv )
{
  CmdLine cmd;

  std::string sInputImage ;
  std::string sOuputFile ;

  int iNbOctave = 4 ;
  int iNbSlicePerOctave = 4 ;

  cmd.add( make_option('i', sInputImage , "input") );
  cmd.add( make_option('o', sOuputFile , "output") );
  cmd.add( make_option('p', iNbOctave , "nb-octave" ) ) ;
  cmd.add( make_option('q', iNbSlicePerOctave , "nb-slice" ) ) ;

  if( argc == 1 )
  {
    std::cerr << "Error : No option given" << std::endl ;
    usage( argv[0] ) ;
  }

  cmd.process(argc, argv);


  if( sInputImage.empty() )
  {
    std::cerr << "Error : input file name empty" << std::endl ;
    usage( argv[0] ) ;
  }

  if( sOuputFile.empty() )
  {
    std::cerr << "Error : output file name empty" << std::endl ;
    usage( argv[0] ) ;
  }

  // Compute base output filename
  const std::string outputBaseName = stlplus::basename_part( sOuputFile ) ;

  Image<unsigned char> src ;
  ReadImage( sInputImage.c_str() , &src ) ;

  Timer t;
  t.reset();

  AKAZEConfig options;
  options.fDesc_factor = 10.f * sqrt(2.f) ;
  AKAZE akaze(src, options);
  akaze.Compute_AKAZEScaleSpace();
  std::vector<AKAZEKeypoint> kpts;
  kpts.reserve(5000);
  akaze.Feature_Detection(kpts);
  
  std::cout << "in "
    << t.elapsedMs() << " msec." << std::endl 
    << t.elapsed() << " sec." << std::endl;

  akaze.Do_Subpixel_Refinement(kpts);

  for (size_t i = 0; i < kpts.size(); ++i)
  {
    AKAZEKeypoint & pt = kpts[i];
    akaze.Compute_Main_Orientation(pt,
      akaze.getSlices()[pt.class_id].Lx,
      akaze.getSlices()[pt.class_id].Ly);
  }
  std::cout << "Found " << kpts.size() << " keypoints" << std::endl;

  for (size_t i = 0; i < kpts.size(); ++i)
  {
    const AKAZEKeypoint & kp = kpts[i];
    float ratio = pow(2.f,kp.octave);
    DrawCircle(kp.x, kp.y, kp.size*2.5, 255, &src);
  }

  WriteImage( (outputBaseName + std::string("_feat.png")).c_str(), src);

  svgDrawer svgStream( src.Width(), src.Height());
  svgStream.drawImage(sInputImage, src.Width(), src.Height());

  //-- Draw features
  for (size_t i=0; i< kpts.size(); ++i)  {
    const AKAZEKeypoint & kp = kpts[i];
    float ratio = pow(2.f,kp.octave);
    svgStream.drawCircle(kp.x, kp.y, kp.size*2.50,
        svgStyle().stroke("yellow", 1.0));

    svgStream.drawLine(
      kp.x, kp.y,
      kp.x + cos(kp.angle)*kp.size*2.5, kp.y + sin(kp.angle)*kp.size*2.5,
      svgStyle().stroke("blue", 1.0));
  }

  // Write the SVG file
  std::ofstream svgFile( (outputBaseName + std::string("_feat.svg")).c_str() );
  svgFile << svgStream.closeSvgFile().str();
  svgFile.close();

  return EXIT_SUCCESS ;
}
