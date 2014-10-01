
// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/image/image.hpp"
#include "openMVG/image_features/image_features_fast_dipole.hpp"

#include "software/VO/CGlWindow.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "openMVG/system/timer.hpp"

#include <stdlib.h>
#include <iostream>

int main(int argc, char **argv)
{
  using namespace std;
  std::cout << "VISUAL ODOMETRY" << std::endl;

  CmdLine cmd;

  std::string sImaDirectory = "";
  cmd.add( make_option('i', sImaDirectory, "imadir") );

  try {
    if (argc == 1) throw std::string("Invalid command line parameter.");
    cmd.process(argc, argv);
  } catch(const std::string& s) {
    std::cerr << "Usage: " << argv[0] << '\n'
    << "[-i|--imadir path] \n"
    << std::endl;

    std::cerr << s << std::endl;
    return EXIT_FAILURE;
  }

  if (sImaDirectory.empty())  {
    std::cerr << "\nIt is an invalid output directory" << std::endl;
    return EXIT_FAILURE;
  }

  if ( !glfwInit() )  {
    return 0;
  }

  //--

  CGlWindow window;
  GLuint text2D;

  //---------------------------------------
  // VO: Visual Odometry
  //---------------------------------------

  // Loop over the image sequence
  //  . detect keypoints
  //  . show them in a openGL window
  //  . track features
  //  . perform VO

  std::vector<std::string> vec_image = stlplus::folder_files( sImaDirectory );
  std::sort(vec_image.begin(), vec_image.end());

  Image<unsigned char> currentImage;

  for (std::vector<std::string>::const_iterator iterFile = vec_image.begin();
    iterFile != vec_image.end(); ++iterFile)
  {
    const std::string sImageFilename = stlplus::create_filespec( sImaDirectory, *iterFile );
    if (openMVG::ReadImage( sImageFilename.c_str(), &currentImage))
    {
      if (window._height < 0)
      {
        // no window created yet, initialize it with the first frame
        window.Init(640, 480, "VisualOdometry--TrackingViewer");
        glGenTextures(1,&text2D);             //allocate the memory for texture
        glBindTexture(GL_TEXTURE_2D,text2D);  //Binding the texture
        glEnable(GL_TEXTURE_2D);              //Enable texture
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_LUMINANCE,
          currentImage.Width(), currentImage.Height(), 0,
          GL_LUMINANCE, GL_UNSIGNED_BYTE, currentImage.data());
      }
      //-- Update the openGL texture
      glTexSubImage2D(GL_TEXTURE_2D, 0,0,0,
        currentImage.Width(), currentImage.Height(),
        GL_LUMINANCE, GL_UNSIGNED_BYTE,
        currentImage.data());

      //-- Draw the current image
      window.SetOrtho(currentImage.Width(), currentImage.Height());
      window.DrawFullScreenTexQuad(currentImage.Width(), currentImage.Height());

      // Clear the depth buffer so the drawn image become the background
	    glClear(GL_DEPTH_BUFFER_BIT);

      //--
	    //-- Feature detection - description
      //--

      using namespace openMVG;
      using namespace openMVG::image_features;
	    image_features::Regions * img_kps = NULL;
	    image_features::Fast_dipole_img_features kp_extractor;
	    kp_extractor.detect_and_describe(currentImage, img_kps);

      const std::vector<PointFeature> kpts = img_kps->getRegionsPositions();
	    std::cout << "Found #" << kpts.size() << std::endl;
	    //--
	    // Draw keypoints
	    glDisable(GL_DEPTH);
	    glPointSize(5.0f);
	    glColor3f(0.f, 0.f, 1.f);
	    glBegin(GL_POINTS);
	    for (size_t i = 0; i < kpts.size(); ++i)
        glVertex2f(kpts[i].x(), kpts[i].y());
      glEnd();

	    //--
	    //-- Feature tracking (todo)
      //--

      if (img_kps) delete img_kps;

      window.Swap();
    }
  }

  // Let last image on the screen
  while(!glfwWindowShouldClose(window._window))
    window.Swap();

  glfwTerminate();
  return 0;
}
