
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Implementation corresponding to the ACCV12 demo:
// Adaptive model estimation, a real time demonstration.
// Pierre Moulon, Pascal Monasse and Renaud Marlet.
// In 11th Asian Confence on Computer Vision (ACCV 2012)

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>

#include <cv.h>
#include <highgui.h>
#include <cvaux.h>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/features2d/features2d.hpp>
#include <opencv2/nonfree/features2d.hpp>
#include <opencv2/legacy/legacy.hpp>

#include "openMVG/numeric/numeric.h"
#include "third_party/vectorGraphics/svgDrawer.hpp"
#include "openMVG/multiview/solver_homography_kernel.hpp"
#include "openMVG/multiview/conditioning.hpp"

using namespace openMVG;

#include "openMVG/robust_estimation/robust_estimator_MaxConsensus.hpp"
#include "openMVG/robust_estimation/score_evaluator.hpp"

#include "openMVG/robust_estimation/robust_estimator_ACRansac.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansacKernelAdaptator.hpp"

using namespace openMVG::robust;


#include <boost/circular_buffer.hpp>
#include <numeric>

using namespace cv;

// Number of the frames processed per second in the application
int fps;

// Estimates the fps of the application
void fpsCalculation(void)
{
  static int64 currentTime, lastTime = cvGetTickCount();
  static int fpsCounter = 0;
  currentTime = cvGetTickCount();
  ++fpsCounter;

  // If 1 second has passed since the last FPS estimation, update the fps
  if (currentTime - lastTime > 1e6 * cvGetTickFrequency()) {
    fps = fpsCounter;
    lastTime = currentTime;
    fpsCounter = 0;
  }
}

int fontFace = CV_FONT_HERSHEY_PLAIN;
double fontScale = 1;
int thickness = 1;

int main(int, char**)
{
  cout << "Press:" << endl
       << " 'q' to quit." << endl;

  cv::VideoCapture cap(1); // open the Left camera
  int w = 640, h = 480;
  cap.set(CV_CAP_PROP_FRAME_WIDTH, w);
  cap.set(CV_CAP_PROP_FRAME_HEIGHT, h);

  if(!cap.isOpened() )  // check if we succeeded
  {
      std::cerr << "ERROR: Could not open cameras." << std::endl;
      return -1;
  }

  std::vector< cv::KeyPoint > kptsStream;
  kptsStream.reserve(5000);
  cv::Mat descriptorsStream;

  bool bDoRansac = true;

  Ptr<FeatureDetector> detector;
  detector = new SiftFeatureDetector;

  Ptr<DescriptorExtractor> extractor;
  extractor = new FREAK(true, true);

  std::vector< cv::KeyPoint > kpts;
  cv::Mat desc;
  std::string filename = string(THIS_SOURCE_DIR) + string("/ACCVPattern.jpg");
  cv::Mat imgA = imread( filename.c_str(), CV_LOAD_IMAGE_GRAYSCALE );
  cv::Mat imgAColor = imread(filename, CV_LOAD_IMAGE_COLOR );
  //Build template reference :
  {
    detector->detect(imgA, kpts);
    extractor->compute(imgA, kpts, desc);
    std::cout << "-------------------\n"
    << "template : " << kpts.size() << " points"
    << "-------------------\n"
    << std::endl;
  }
  cv::Mat sidebyside;
  sidebyside.create(std::max(imgA.rows, h)+10, imgA.cols+w+10, CV_8UC3);

  //copy template
  drawKeypoints(imgAColor, kpts, imgAColor);
  imgAColor.copyTo(sidebyside.colRange(0,0+imgAColor.cols).rowRange(0,0+imgAColor.rows));

  BFMatcher matcher(NORM_HAMMING, true);
  std::vector<cv::DMatch> vec_matches;

  //--------
  //-- Statistics
  boost::circular_buffer<double> cb_ACRANSAC_Threshold(100,0);
  boost::circular_buffer<int> cb_ACRANSAC_found(100,0);
  boost::circular_buffer<int> cb_RANSAC_found(100,0);

  boost::circular_buffer<double> cb_ACRANSAC_inlierPourcent(100,0);
  boost::circular_buffer<double> cb_RANSAC_inlierPourcent(100,0);

  cv::namedWindow("WEBCAM STREAM",CV_WINDOW_AUTOSIZE);
  for(;;)
  {
      bool isValid = true;

      fpsCalculation();

      cv::Mat frame, grayFrame;

      try
      {
        cap >> frame; // get a new frame from left camera
        cv::cvtColor(frame, grayFrame, CV_RGB2GRAY);
        detector->detect(grayFrame, kptsStream);
        extractor->compute( grayFrame, kptsStream, descriptorsStream );
        //std::cout << "KpFound : \t" << kptsStream.size() << std::endl;
      }
      catch( cv::Exception& e )
      {
        std::cout << "An exception occurred. Ignoring frame. " << e.err << std::endl;
        isValid = false;
      }

      if (isValid)
      {
        try
        {
          static char text[256];
          sprintf(text, "FPS: %d", fps);

          int baseline = 0;
          Size textSize = getTextSize(text, fontFace,
                            fontScale, thickness, &baseline);
          putText(frame, text, cv::Point(10, 30), fontFace, fontScale,
            Scalar::all(255), thickness, 8);

          drawKeypoints(frame, kptsStream, frame);

          matcher.match(desc, descriptorsStream, vec_matches);
          //std::cout << "found : " << vec_matches.size() << "matches" << std::endl;

          frame.copyTo(sidebyside.colRange(imgA.cols,imgA.cols+frame.cols).rowRange(0,frame.rows));

          //-------
          //-- Robust Model estimation :
          //-------
          openMVG::Mat xA(2,vec_matches.size()), xB(2, vec_matches.size());
          const std::vector<cv::KeyPoint >& keypointsA = kpts;
          const std::vector<cv::KeyPoint >& keypointsB = kptsStream;
          for (size_t i = 0; i < vec_matches.size(); ++i)
          {
            const DMatch& match = vec_matches[i];
            xA.col(i) = Vec2(keypointsA[match.queryIdx].pt.x, keypointsA[match.queryIdx].pt.y);
            xB.col(i) = Vec2(keypointsB[match.trainIdx].pt.x, keypointsB[match.trainIdx].pt.y);
          }

          //-- Homography robust estimation
          std::vector<size_t> vec_inliers, vec_inliersRansac;
          Mat3 H;
          Mat3 Hransac;
          double thresholdransac = 2.0;
          double thresholdH = 4.0;
          double NFAH;
          {
            typedef ACKernelAdaptor<
              openMVG::homography::kernel::FourPointSolver,
              openMVG::homography::kernel::AsymmetricError,
              UnnormalizerI,
              Mat3>
            KernelType;

            KernelType kernel(
                  xA, imgA.cols, imgA.rows,
                  xB, frame.cols, frame.rows,
                  false); // configure as point to point error model.

            double error = std::numeric_limits<double>::infinity();
            const std::pair<double,double> ACRansacOut =
              ACRANSAC(kernel, vec_inliers, 1024, &H, error, false);
            thresholdH = ACRansacOut.first;
            NFAH = ACRansacOut.second;
          }
          {
            typedef homography::kernel::Kernel KernelType;
            KernelType kernel(xA, xB);
            if (bDoRansac)
              Hransac = MaxConsensus(kernel, ScorerEvaluator<KernelType>(Square(thresholdransac)), &vec_inliersRansac, 1024);
            if (vec_inliersRansac.size()<12)
              thresholdransac = 900;
          }

          std::vector<char> vec_matchesMask(vec_matches.size(), 0);
          for (size_t i = 0; i < vec_inliers.size(); ++i)
            vec_matchesMask[ vec_inliers[i] ] = 1;

          static const int Lthickness = 3;

          Scalar ACColor(255,255,255);
          Scalar RColor(255,0,255);
          bool bDrawRansac = thresholdransac < 40;
          bool bDrawACRansac = thresholdH < 15 && vec_inliers.size() > 12;

          if (bDrawACRansac)
            drawMatches(imgA, kpts, frame, kptsStream, vec_matches, sidebyside,
                      Scalar(0,255,0), Scalar(0,255,255), vec_matchesMask);

          if (bDrawRansac)
          {
            cb_RANSAC_found.push_back(1);
            cb_RANSAC_inlierPourcent.push_back(vec_inliersRansac.size()/float(vec_matches.size()));
            //Draw warp
            Vec2 x0(0,0), x1(imgA.cols, 0), x2(imgA.cols, imgA.rows), x3(0, imgA.rows);
            Vec3 x00 = Hransac * Vec3(x0(0), x0(1), 1);
            Vec3 x11 = Hransac * Vec3(x1(0), x1(1), 1);
            Vec3 x22 = Hransac * Vec3(x2(0), x2(1), 1);
            Vec3 x33 = Hransac * Vec3(x3(0), x3(1), 1);
            x0 = x00.head<2>() / x00(2);
            x1 = x11.head<2>() / x11(2);
            x2 = x22.head<2>() / x22(2);
            x3 = x33.head<2>() / x33(2);

            line(sidebyside, Point(x0(0)+imgA.cols, x0(1)), Point(x1(0)+imgA.cols, x1(1)), RColor, Lthickness);
            line(sidebyside, Point(x1(0)+imgA.cols, x1(1)), Point(x2(0)+imgA.cols, x2(1)), RColor, Lthickness);
            line(sidebyside, Point(x2(0)+imgA.cols, x2(1)), Point(x3(0)+imgA.cols, x3(1)), RColor, Lthickness);
            line(sidebyside, Point(x3(0)+imgA.cols, x3(1)), Point(x0(0)+imgA.cols, x0(1)), RColor, Lthickness);
          }
          else{
            cb_RANSAC_found.push_back(0);
            cb_RANSAC_inlierPourcent.push_back(0);
          }

          if (bDrawACRansac)
          {
            cb_ACRANSAC_Threshold.push_back(thresholdH);
            cb_ACRANSAC_found.push_back(1);
            cb_ACRANSAC_inlierPourcent.push_back(vec_inliers.size()/float(vec_matches.size()));

            //Draw warp
            Vec2 x0(0,0), x1(imgA.cols, 0), x2(imgA.cols, imgA.rows), x3(0, imgA.rows);
            Vec3 x00 = H * Vec3(x0(0), x0(1), 1);
            Vec3 x11 = H * Vec3(x1(0), x1(1), 1);
            Vec3 x22 = H * Vec3(x2(0), x2(1), 1);
            Vec3 x33 = H * Vec3(x3(0), x3(1), 1);
            x0 = x00.head<2>() / x00(2);
            x1 = x11.head<2>() / x11(2);
            x2 = x22.head<2>() / x22(2);
            x3 = x33.head<2>() / x33(2);

            line(sidebyside, Point(x0(0)+imgA.cols, x0(1)), Point(x1(0)+imgA.cols, x1(1)), ACColor, Lthickness);
            line(sidebyside, Point(x1(0)+imgA.cols, x1(1)), Point(x2(0)+imgA.cols, x2(1)), ACColor, Lthickness);
            line(sidebyside, Point(x2(0)+imgA.cols, x2(1)), Point(x3(0)+imgA.cols, x3(1)), ACColor, Lthickness);
            line(sidebyside, Point(x3(0)+imgA.cols, x3(1)), Point(x0(0)+imgA.cols, x0(1)), ACColor, Lthickness);

            std::stringstream stext;
            stext << " inlier/outlier ratio: " << vec_inliers.size()/float(vec_matches.size());

            baseline = 0;
            textSize = getTextSize(stext.str().c_str(), fontFace,
                              fontScale, thickness, &baseline);
            putText(sidebyside, stext.str().c_str(), cv::Point(imgA.cols+10, 60), fontFace, fontScale,
              Scalar::all(255), thickness, 8);
          }
          else
          {
            thresholdH = 0;
            vec_inliers.resize(0);
            cb_ACRANSAC_Threshold.push_back(0.0);
            cb_ACRANSAC_found.push_back(0);
            cb_ACRANSAC_inlierPourcent.push_back(0);
          }


          cv::imshow("WEBCAM STREAM", sidebyside);

          //-- Display stats:

          cv::Mat plotImg = cv::Mat::zeros( 100 , cb_ACRANSAC_Threshold.size(), CV_8UC3);
          cv::Point * arrayA = new cv::Point[100];
          for (int i=0; i < cb_ACRANSAC_Threshold.size(); ++i){
            arrayA[i] = Point(i, 100-cb_ACRANSAC_Threshold[i]*5.0);
          }
          int nCurvePts[1]={cb_ACRANSAC_Threshold.size()};
          cv::polylines(plotImg, (const cv::Point**)&arrayA, nCurvePts, 1, false, ACColor, 1);

          std::stringstream stext;
          stext << "T = :" << thresholdH;
          baseline = 0;
          textSize = getTextSize(stext.str().c_str(), fontFace, fontScale, thickness, &baseline);
          putText(plotImg, stext.str().c_str(), cv::Point(10, 100-80), fontFace, fontScale, Scalar::all(255), thickness, 8);

          cv::imshow("AC-RANSAC Threshold", plotImg);

          //Display found/not found
          plotImg = cv::Mat::zeros( 100 , cb_ACRANSAC_Threshold.size(), CV_8UC3);
          for (int i=0; i < cb_ACRANSAC_Threshold.size(); ++i){
            arrayA[i] = Point(i, cb_ACRANSAC_found[i]*8.0);
          }
          cv::polylines(plotImg, (const cv::Point**)&arrayA, nCurvePts, 1, false, ACColor, 2);

          for (int i=0; i < cb_RANSAC_found.size(); ++i){
            arrayA[i] = Point(i, cb_RANSAC_found[i]*8.0 + 50);
          }
          cv::polylines(plotImg, (const cv::Point**)&arrayA, nCurvePts, 1, false, RColor, 2);

          stext.str("");
          stext << "#found :" << std::accumulate(cb_ACRANSAC_found.begin(), cb_ACRANSAC_found.end(), 0);

          baseline = 0;
          textSize = getTextSize(stext.str().c_str(), fontFace,
                            fontScale, thickness, &baseline);
          putText(plotImg, stext.str().c_str(), cv::Point(1, 30), fontFace, fontScale,
            Scalar::all(255), thickness, 8);

          stext.str("");
          stext << "#found :" << std::accumulate(cb_RANSAC_found.begin(), cb_RANSAC_found.end(), 0);
          baseline = 0;
          textSize = getTextSize(stext.str().c_str(), fontFace,
                            fontScale, thickness, &baseline);
          putText(plotImg, stext.str().c_str(), cv::Point(1, 80), fontFace, fontScale,
            Scalar::all(255), thickness, 8);

          cv::imshow("Found/NotFound", plotImg);

          //Display inliers pourcentage
          plotImg = cv::Mat::zeros( 104 , cb_ACRANSAC_inlierPourcent.size(), CV_8UC3);
          for (int i=0; i < cb_ACRANSAC_inlierPourcent.size(); ++i){
            arrayA[i] = Point(i, 100-cb_ACRANSAC_inlierPourcent[i]*100.0);
          }
          cv::polylines(plotImg, (const cv::Point**)&arrayA, nCurvePts, 1, false, ACColor, 1);

          for (int i=0; i < cb_RANSAC_inlierPourcent.size(); ++i){
            arrayA[i] = Point(i, 100-cb_RANSAC_inlierPourcent[i]*100.0);
          }
          cv::polylines(plotImg, (const cv::Point**)&arrayA, nCurvePts, 1, false, RColor, 1);

          stext.str("");
          stext << "AC Inlier";
          baseline = 0;
          textSize = getTextSize(stext.str().c_str(), fontFace,
                            fontScale, thickness, &baseline);
          putText(plotImg, stext.str().c_str(), cv::Point(1, 20), fontFace, fontScale,
            Scalar::all(255), thickness, 8);
          stext.str("");
          stext << (int)(vec_inliers.size()/float(vec_matches.size())*100) << "%";
          baseline = 0;
          textSize = getTextSize(stext.str().c_str(), fontFace,
                            fontScale, thickness, &baseline);
          putText(plotImg, stext.str().c_str(), cv::Point(1, 40), fontFace, fontScale,
            Scalar::all(255), thickness, 8);

          cv::imshow("InlierPourcent", plotImg);

          delete [] arrayA;

        }
        catch( cv::Exception& e )
        {
          std::cout << "An exception occurred. Ignoring frame. " << e.err << std::endl;
        }
      }

      bool bQuit = false;
      const char key = cvWaitKey(30);

      switch (key) {
      case 'q': case 'Q': bQuit = true; break;
      case 'c': case 'C':
        {
          fill_n(cb_ACRANSAC_Threshold.begin(),cb_ACRANSAC_Threshold.size(),0);
          fill_n(cb_ACRANSAC_found.begin(),100,0);
          fill_n(cb_RANSAC_found.begin(),100,0);
          fill_n(cb_ACRANSAC_inlierPourcent.begin(),100,0.0);
          fill_n(cb_RANSAC_inlierPourcent.begin(),100,0.0);
        }
        break;
      case 'r': case 'R': bDoRansac = !bDoRansac; break;
      }
      if (bQuit)
        break;
  }
  // the camera object will be automatically released in VideoCapture destructor
  return 0;
}

