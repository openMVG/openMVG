
// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef TRACKER_VO_HPP
#define TRACKER_VO_HPP

// OpenCV Includes
#include <opencv2/core/eigen.hpp> //To Convert Eigen matrix to cv matrix
#include <opencv2/features2d/features2d.hpp>
#include <opencv2/video/tracking.hpp>

#include <vector>

namespace openMVG  {
namespace VO  {

/// OpenCV based KLT tracker
struct Tracker_opencv_KLT
{
  // data for tracking
  cv::Mat _current_img, _prev_img;
  std::vector<cv::Point2f> _prevPts, _nextPts;
  std::vector<float> _error;

  /// Try to track current point set in the provided image
  /// return false when tracking failed (=> to send frame to relocalization)
  bool track(
    const Image<unsigned char> & ima,
    const std::vector<PointFeature> & pt_to_track,
    std::vector<PointFeature> & pt_tracked,
    std::vector<unsigned char> & status)
  {
    cv::eigen2cv(ima.GetMat(), _current_img);
    if (!pt_to_track.empty())
    {
      _prevPts.resize(pt_to_track.size());
      _nextPts.resize(pt_to_track.size());

      for (size_t i=0; i < pt_to_track.size(); ++i)
      {
        _prevPts[i].x = pt_to_track[i].x();
        _prevPts[i].y = pt_to_track[i].y();
      }

      cv::calcOpticalFlowPyrLK(_prev_img, _current_img, _prevPts, _nextPts, status, _error);

      for (size_t i=0; i < _nextPts.size(); ++i)
      {
        pt_tracked[i].coords() = Vec2f(_nextPts[i].x, _nextPts[i].y);
      }
    }
    // swap frame for next tracking iteration
    _current_img.copyTo(_prev_img);
    return false;
  }

  // suggest new count feature point for tracking
  bool detect(
    const Image<unsigned char> & ima,
    std::vector<PointFeature> & pt_to_track,
    const size_t count)
  {
    cv::Mat current_img;
    cv::eigen2cv(ima.GetMat(), current_img);
    std::vector<cv::KeyPoint> m_nextKeypoints;
    cv::Ptr<cv::FeatureDetector> m_detector = cv::FeatureDetector::create("GridFAST");
    m_detector->detect(current_img, m_nextKeypoints);
    if (m_nextKeypoints.size() > count)
    {
      // shuffle to avoid to sample only in one bucket
      std::random_shuffle(m_nextKeypoints.begin(), m_nextKeypoints.end());
      m_nextKeypoints.resize(count);

      pt_to_track.resize(count);
      for (size_t i = 0; i  < count; ++i)
        pt_to_track[i] = (PointFeature(m_nextKeypoints[i].pt.x, m_nextKeypoints[i].pt.y));

      return true;
    }
    else
      return false;
  }
};

} // namespace VO
} // namespace openMVG


#endif // TRACKER_VO_HPP
