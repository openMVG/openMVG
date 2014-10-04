
// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/image/image.hpp"
#include "openMVG/features/feature.hpp"

#include "software/VO/CGlWindow.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "openMVG/system/timer.hpp"

// OpenCV Includes
#include "opencv2/core/eigen.hpp" //To Convert Eigen matrix to cv matrix
#include <cv.h>

#include <stdlib.h>
#include <deque>
#include <iostream>

using namespace openMVG;

// Struct to store a partial PointFeature trajectory
struct STrajectory
{
  STrajectory() : _bAlive(false) {}
  STrajectory(const PointFeature & pos) : _bAlive(true) { _history.push_back(pos); }

  std::deque<PointFeature> _history;
  bool _bAlive;
};

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

  if (sImaDirectory.empty() || !stlplus::is_folder(sImaDirectory))  {
    std::cerr << "\nIt is an invalid input directory" << std::endl;
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

  // OpenCV Data (required for KLT tracking):
  cv::Mat current_img, prev_img;
  std::vector<cv::Point2f> _prevPts, _nextPts;
  cv::Ptr<cv::FeatureDetector> m_detector = cv::FeatureDetector::create("GridFAST");
  const size_t cardinalFeature = 500;
  std::vector<STrajectory> vec_trajectory(cardinalFeature);

  size_t iterCount = 0;
  for (std::vector<std::string>::const_iterator iterFile = vec_image.begin();
    iterFile != vec_image.end(); ++iterFile, ++iterCount)
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

      //Convert image to OpenCV data

      cv::eigen2cv(currentImage.GetMat(), current_img);

      std::vector<unsigned char> status;
      std::vector<float> error;

      if (_prevPts.size() > 0)  {
        cv::calcOpticalFlowPyrLK(prev_img, current_img, _prevPts, _nextPts, status, error);
      }

      // Count the number of tracked features
      const size_t countTracked = accumulate(status.begin(), status.end(), 0);

      // Update the list of tracked features (the not tracked one will be updated with new feature later)
      std::vector<cv::Point2f> trackedPts(cardinalFeature);
      for (size_t i=0; i<status.size(); ++i)  {
        if (status[i])  trackedPts[i] = _nextPts[i];
      }

      // Update trajectories
      for(unsigned i = 0; i < status.size(); ++i) {
        if (status[i])  {
          // if tracked
          const PointFeature a(_prevPts[i].x, _prevPts[i].y);
          const PointFeature b(_nextPts[i].x, _nextPts[i].y);
          //std::cout << DistanceL2(a.coords(), b.coords()) << " ";

          vec_trajectory[i]._history.push_back(b);
          vec_trajectory[i]._bAlive = true;
          if (vec_trajectory[i]._history.size() > 10) vec_trajectory[i]._history.pop_front();
        }
        else  {
          // feature was not re-detected
          vec_trajectory[i]._history.clear();
          vec_trajectory[i]._bAlive = false;
        }
      }

      //--
      // Draw feature trajectories
      //--
      glColor3f(0.f, 1.f, 0.f);
      glLineWidth(2.f);
      for(size_t i = 0; i < vec_trajectory.size(); ++i)
      {
        if (vec_trajectory[i]._bAlive && vec_trajectory[i]._history.size() > 1)
        {
          glBegin(GL_LINE_STRIP);
          for (unsigned j=0; j<vec_trajectory[i]._history.size(); j++)
          {
            const PointFeature & p0 = vec_trajectory[i]._history[j];
            glVertex2f(p0.x(), p0.y());
          }
          glEnd();
        }
      }
      glFlush();

      // Do we have to add new feature (if some track have been cut)
      const bool bTooFewFeature = countTracked < cardinalFeature;
      if (bTooFewFeature)
      {
        std::vector<cv::KeyPoint> m_nextKeypoints;
        m_detector->detect(current_img, m_nextKeypoints);
        const size_t pointsToDetect = cardinalFeature - countTracked;
        if (m_nextKeypoints.size() > pointsToDetect)
        {
          // shuffle to avoid to sample only in one bucket
          std::random_shuffle(m_nextKeypoints.begin(), m_nextKeypoints.end());
          m_nextKeypoints.resize(pointsToDetect);
        }
        std::cout << "#features added: " << m_nextKeypoints.size() << std::endl;
        size_t j = 0;
        for(size_t i = 0; i < vec_trajectory.size(); ++i)
        {
          if (!vec_trajectory[i]._bAlive) // update the feature
          {
            trackedPts[i] = m_nextKeypoints[j].pt;
            vec_trajectory[i] = STrajectory(PointFeature(m_nextKeypoints[j].pt.x, m_nextKeypoints[j].pt.y));
            ++j;
          }
        }
      }
      // Swap image and tracked point for future iteration
      _prevPts.swap(trackedPts);
      current_img.copyTo(prev_img);

      window.Swap(); // Swap openGL buffer
    }
  }

  // Let last image on the screen
  while(!glfwWindowShouldClose(window._window))
    window.Swap();

  glfwTerminate();
  return 0;
}
