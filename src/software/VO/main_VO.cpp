// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/image/image_io.hpp"
#include "openMVG/features/feature.hpp"

#include "software/VO/CGlWindow.hpp"
#include "software/VO/Monocular_VO.hpp"
#include "software/VO/Tracker.hpp"
#if defined HAVE_OPENCV
#include "software/VO/Tracker_opencv_klt.hpp"
#endif

#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <cstdlib>
#include <iostream>

using namespace openMVG;

int main(int argc, char **argv)
{
  using namespace std;
  std::cout << "VISUAL ODOMETRY -- Tracking demo --" << std::endl;

  CmdLine cmd;

  std::string sImaDirectory = "";
  std::string sOutFile = "";
  unsigned int uTracker = 0;
  unsigned int uTrackerPointCount = 1500;

  cmd.add( make_option('i', sImaDirectory, "imadir") );
  cmd.add( make_option('t', uTracker, "tracker") );
  cmd.add( make_option('o', sOutFile, "output_file") );
  cmd.add( make_option('p', uTrackerPointCount, "point_count") );
  cmd.add( make_switch('d', "disable_tracking_display") );

  try {
    if (argc == 1) throw std::string("Invalid command line parameter.");
    cmd.process(argc, argv);
  } catch (const std::string& s) {
    std::cerr << "Usage: " << argv[0] << '\n'
    << "[-i|--imadir path] \n"
    << "[-o|--output_file path] The output tracking saved as a sfm_data file (landmark observations).\n"
    << "[-t|--tracker Used tracking interface] \n"
    << "\t 0: Feature matching based tracking (default); Fast detector + Dipole descriptor, \n"
#if defined HAVE_OPENCV
    << "\t 1: Feature tracking based tracking; Fast + KLT pyramidal tracking. \n"
#endif
    << "[-p|--point_count] Number of points to track. (default: " << uTrackerPointCount << ")\n"
    << "[-d|--disable_tracking_display] Disable tracking display \n"
    << std::endl;

    std::cerr << s << std::endl;
    return EXIT_FAILURE;
  }

   std::cout << " You called : " << std::endl
            << argv[0] << std::endl
            << "--imageDirectory " << sImaDirectory << std::endl
            << "--output_file " << sOutFile << std::endl
            << "--point_count " << uTrackerPointCount << std::endl
            << "--tracker " << uTracker << std::endl
            << "--disable_tracking_display " << static_cast<int>(cmd.used('d')) << std::endl;

  if (sImaDirectory.empty() || !stlplus::is_folder(sImaDirectory))
  {
    std::cerr << "\nIt is an invalid input directory" << std::endl;
    return EXIT_FAILURE;
  }

  if (sOutFile.empty())
  {
    std::cerr << "\nPlease use a valid filename for the output_file option (i.e: <PATH>/sfm_data.bin)." << std::endl;
    return EXIT_FAILURE;
  }

  if ( !glfwInit() )
  {
    return EXIT_FAILURE;
  }

  const bool disable_tracking_display = cmd.used('d');

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

  std::vector<std::string> vec_image = stlplus::folder_files(sImaDirectory);
  // clean invalid image file
  {
    std::vector<std::string> vec_image_;
    for (size_t i = 0; i < vec_image.size(); ++i)
    {
      if (openMVG::image::GetFormat(vec_image[i].c_str()) != openMVG::image::Unknown)
        vec_image_.push_back(vec_image[i]);
    }
    vec_image_.swap(vec_image);
  }
  std::sort(vec_image.begin(), vec_image.end());

  image::Image<unsigned char> currentImage;

  using namespace openMVG::VO;

  // Initialize the tracker interface
  std::unique_ptr<Abstract_Tracker> tracker_ptr;
  switch (uTracker)
  {
    case 0:
      tracker_ptr.reset(new Tracker_fast_dipole);
    break;
#if defined HAVE_OPENCV
    case 1:
      tracker_ptr.reset(new Tracker_opencv_KLT);
    break;
#endif
    default:
    std::cerr << "Unknow tracking method" << std::endl;
    return EXIT_FAILURE;
  }
  if (!tracker_ptr)
  {
    std::cerr << "Cannot instantiate the Feature tracking interface" << std::endl;
    return EXIT_FAILURE;
  }

  // Initialize the monocular tracking framework
  VO_Monocular monocular_vo(tracker_ptr.get(), uTrackerPointCount);

  size_t frameId = 0;
  for (std::vector<std::string>::const_iterator iterFile = vec_image.begin();
    iterFile != vec_image.end(); ++iterFile, ++frameId)
  {
    const std::string sImageFilename = stlplus::create_filespec( sImaDirectory, *iterFile );
    if (openMVG::image::ReadImage( sImageFilename.c_str(), &currentImage))
    {
      if (window._height < 0)
      {
        // no window created yet, initialize it with the first frame

        const double aspect_ratio = currentImage.Width() / (double)currentImage.Height();
        window.Init(640, 640 / aspect_ratio, "VisualOdometry--TrackingViewer");
        glGenTextures(1, &text2D);             //allocate the memory for texture
        glBindTexture(GL_TEXTURE_2D, text2D); //Binding the texture
        glEnable(GL_TEXTURE_2D);              //Enable texture
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_LUMINANCE,
          currentImage.Width(), currentImage.Height(), 0,
          GL_LUMINANCE, GL_UNSIGNED_BYTE, currentImage.data());
      }
      glBindTexture(GL_TEXTURE_2D, text2D); //Binding the texture
      glEnable(GL_TEXTURE_2D);              //Enable texture
      //-- Update the openGL texture with the current frame pixel values
      glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0,
        currentImage.Width(), currentImage.Height(),
        GL_LUMINANCE, GL_UNSIGNED_BYTE,
        currentImage.data());

      //-- Draw the current image
      window.SetOrtho(currentImage.Width(), currentImage.Height());
      window.DrawFullScreenTexQuad(currentImage.Width(), currentImage.Height());
      glDisable(GL_TEXTURE_2D);

      // Clear the depth buffer so the drawn image become the background
      glClear(GL_DEPTH_BUFFER_BIT);
      glDisable(GL_LIGHTING);

      //--
      //-- Feature tracking
      //    . track features
      //    . if some tracks are cut, detect and insert new features
      //--
      monocular_vo.nextFrame(currentImage, frameId);

      //--
      // Draw feature trajectories
      //--
      glColor3f(0.f, 1.f, 0.f);
      glLineWidth(2.f);

      if (!disable_tracking_display)
      {
        for (size_t idx = 0; idx < monocular_vo.landmark_.size(); ++idx)
        {
          if (std::find(monocular_vo.trackedLandmarkIds_.cbegin(),
                        monocular_vo.trackedLandmarkIds_.cend(), idx)
              == monocular_vo.trackedLandmarkIds_.cend())
            continue;

          const Landmark & landmark = monocular_vo.landmark_[idx];
          if (landmark.obs_.back().frameId_ == frameId && landmark.obs_.size() > 1 )
          {
            const std::deque<Measurement> & obs = landmark.obs_;

            std::deque<Measurement>::const_reverse_iterator
              iter = obs.rbegin(),
              iterEnd = obs.rend();

            int limit = 10;
            glBegin(GL_LINE_STRIP);
            glColor3f(0.f, 1.f, 0.f);
            for (; iter != iterEnd && limit >=0; ++iter, --limit)
            {
              const Vec2f & p0 = iter->pos_;
              glVertex2f(p0(0), p0(1));
            }
            glEnd();

            // draw the current tracked point
            {
              std::deque<Measurement>::const_reverse_iterator iter = obs.rbegin();
              glPointSize(4.0f);
              glBegin(GL_POINTS);
              glColor3f(1.f, 1.f, 0.f); // Yellow
              const Vec2f & p0 = iter->pos_;
              glVertex2f(p0(0), p0(1));
              glEnd();
            }
          }
          else // Draw the new initialized point
          {
            glPointSize(10.0f);
            if ( landmark.obs_.size() == 1 )
            {
              const std::deque<Measurement> & obs = landmark.obs_;
              glBegin(GL_POINTS);
              glColor3f(0.f, 0.f, 1.f); // Blue
              std::deque<Measurement>::const_iterator iter = obs.begin();
              const Vec2f & p0 = iter->pos_;
              glVertex2f(p0(0), p0(1));
              glEnd();
            }
          }
        }
      }

      glFlush();

      window.Swap(); // Swap openGL buffer
    }
  }

  openMVG::sfm::SfM_Data sfm_data;
  ConvertVOLandmarkToSfMDataLandmark(monocular_vo.landmark_, sfm_data.structure);
  std::cout << "Found SFM #landmarks: " << sfm_data.structure.size() << std::endl;
  if (!Save(sfm_data, sOutFile, openMVG::sfm::ESfM_Data(openMVG::sfm::ALL)))
    return EXIT_FAILURE;

  glfwTerminate();
  return 0;
}
