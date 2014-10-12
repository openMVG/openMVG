
// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/image/image.hpp"
#include "openMVG/features/feature.hpp"

#include "software/VO/CGlWindow.hpp"
#include "software/VO/Monocular_VO.hpp"
#include "software/VO/Tracker.hpp"

#include "openMVG/split/split.hpp"
#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <stdlib.h>
#include <iostream>

using namespace openMVG;

/// Check that Kmatrix is a string like "f;0;ppx;0;f;ppy;0;0;1"
/// With f,ppx,ppy as valid numerical value
bool checkIntrinsicStringValidity(const std::string & Kmatrix, Mat3 & K)
{
  std::vector<std::string> vec_str;
  split( Kmatrix, ";", vec_str );
  if (vec_str.size() != 9)  {
    std::cerr << "\n Missing ';' character" << std::endl;
    return false;
  }
  // Check that all K matrix value are valid numbers
  for (size_t i = 0; i < vec_str.size(); ++i) {
    double readvalue = 0.0;
    std::stringstream ss;
    ss.str(vec_str[i]);
    if (! (ss >> readvalue) )  {
      std::cerr << "\n Used an invalid not a number character" << std::endl;
      return false;
    }
    *(K.data()+i) = readvalue;
  }
  return true;
}

int main(int argc, char **argv)
{
  using namespace std;
  std::cout << "VISUAL ODOMETRY" << std::endl;

  CmdLine cmd;

  std::string
    sImaDirectory = "",
    sKmatrix = "";

  cmd.add( make_option('i', sImaDirectory, "imadir") );
  cmd.add( make_option('k', sKmatrix, "intrinsics") );

  try {
    if (argc == 1) throw std::string("Invalid command line parameter.");
    cmd.process(argc, argv);
  } catch(const std::string& s) {
    std::cerr << "Usage: " << argv[0] << '\n'
    << "[-i|--imadir path] \n"
    << "[-k|--intrinsics] Kmatrix: \"f;0;ppx;0;f;ppy;0;0;1\""
    << std::endl;

    std::cerr << s << std::endl;
    return EXIT_FAILURE;
  }

   std::cout << " You called : " <<std::endl
            << argv[0] << std::endl
            << "--imageDirectory " << sImaDirectory << std::endl
            << "--intrinsics " << sKmatrix << std::endl;

  if (sImaDirectory.empty() || !stlplus::is_folder(sImaDirectory))  {
    std::cerr << "\nIt is an invalid input directory" << std::endl;
    return EXIT_FAILURE;
  }

  Mat3 K;
  if (!checkIntrinsicStringValidity(sKmatrix, K))
  {
    std::cerr << "\nInvalid K matrix input" << std::endl;
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
  //  . show them in an openGL window
  //  . track features
  //  . perform VO

  std::vector<std::string> vec_image = stlplus::folder_files( sImaDirectory );
  std::sort(vec_image.begin(), vec_image.end());

  Image<unsigned char> currentImage;

  using namespace openMVG::VO;
  VO_Monocular<Tracker_opencv_KLT> monocular_vo;

  size_t frameId = 0;
  for (std::vector<std::string>::const_iterator iterFile = vec_image.begin();
    iterFile != vec_image.end(); ++iterFile, ++frameId)
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
	    glDisable(GL_LIGHTING);

      //--
	    //-- Feature tracking (KLT)
	    //    . track features
	    //    . if some track are cut, detect and insert new features
      //--

      monocular_vo.nextFrame(currentImage, frameId);

      //--
      // Draw feature trajectories
      //--
      glColor3f(0.f, 1.f, 0.f);
      glLineWidth(2.f);

      for (size_t idx = 0; idx < monocular_vo._landmark.size(); ++idx)
      {
        if (std::find(monocular_vo._trackedLandmarkIds.begin(), monocular_vo._trackedLandmarkIds.end(), idx) == monocular_vo._trackedLandmarkIds.end())
          continue;

        const Landmark & landmark = monocular_vo._landmark[idx];
        if (landmark._obs.back()._frameId == frameId && landmark._obs.size() > 1 )
        {
          const std::deque<Measurement> & obs = landmark._obs;
          glBegin(GL_LINE_STRIP);
          std::deque<Measurement>::const_reverse_iterator iter = obs.rbegin();
          std::deque<Measurement>::const_reverse_iterator iterEnd = obs.rend();

          int limit = 10;
          for (; iter != iterEnd && limit>=0; ++iter, --limit)
          {
            const Vec2f & p0 = iter->_pos;
            glVertex2f(p0(0), p0(1));
          }
          glEnd();
        }

      }
      glFlush();

      window.Swap(); // Swap openGL buffer
    }
  }

  // Let last image on the screen
  while(!glfwWindowShouldClose(window._window))
    window.Swap();

  glfwTerminate();
  return 0;
}
