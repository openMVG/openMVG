// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef TRACKER_OPENCV_KLT_HPP
#define TRACKER_OPENCV_KLT_HPP

// OpenCV Includes
#include <opencv2/core/eigen.hpp> //To Convert Eigen matrix to cv matrix
#include <opencv2/features2d/features2d.hpp>
#include <opencv2/video/tracking.hpp>

#include <software/VO/Abstract_Tracker.hpp>

#include <vector>

namespace openMVG  {
namespace VO  {

/// OpenCV based KLT tracker
struct Tracker_opencv_KLT : public Abstract_Tracker
{
  // data for tracking
  cv::Mat current_img_, prev_img_;
  std::vector<cv::Point2f> prevPts_, nextPts_;
  std::vector<float> error_;

  /// Try to track current point set in the provided image
  /// return false when tracking failed (=> to send frame to relocalization)
  bool track
  (
    const image::Image<unsigned char> & ima,
    const std::vector<features::PointFeature> & pt_to_track,
    std::vector<features::PointFeature> & pt_tracked,
    std::vector<bool> & status
  ) override
  {
    cv::eigen2cv(ima.GetMat(), current_img_);
    if (!pt_to_track.empty())
    {
      prevPts_.resize(pt_to_track.size());
      nextPts_.resize(pt_to_track.size());

      for (size_t i=0; i < pt_to_track.size(); ++i)
      {
        prevPts_[i].x = pt_to_track[i].x();
        prevPts_[i].y = pt_to_track[i].y();
      }

      std::vector<unsigned char> status_uchar;
      cv::calcOpticalFlowPyrLK(prev_img_, current_img_, prevPts_, nextPts_, status_uchar, error_);
      status.assign(status_uchar.begin(), status_uchar.end());

      for (size_t i=0; i < nextPts_.size(); ++i)
      {
        pt_tracked[i].coords() << nextPts_[i].x, nextPts_[i].y;
      }
    }
    // swap frame for next tracking iteration
    current_img_.copyTo(prev_img_);

    const size_t tracked_point_count = std::accumulate(status.begin(), status.end(), 0);
    return (tracked_point_count != 0);
  }

  // suggest new feature point for tracking (count point are kept)
  bool detect
  (
    const image::Image<unsigned char> & ima,
    std::vector<features::PointFeature> & pt_to_track,
    const size_t count
  ) const override
  {
    cv::Mat current_img;
    cv::eigen2cv(ima.GetMat(), current_img);
    std::vector<cv::KeyPoint> m_nextKeypoints;

    cv::Ptr<cv::FeatureDetector> m_detector = cv::GFTTDetector::create(count);
    if (m_detector == NULL)
      return false;

    m_detector->detect(current_img, m_nextKeypoints);

    if (m_nextKeypoints.size() >= count)
    {
      // shuffle to avoid to sample only in one bucket
      std::mt19937 gen(std::mt19937::default_seed);
      std::shuffle(m_nextKeypoints.begin(), m_nextKeypoints.end(), gen);
    }
    const size_t kept_kp_count =  std::min(m_nextKeypoints.size(), count);
    m_nextKeypoints.resize(kept_kp_count);

    pt_to_track.resize(kept_kp_count);
    for (size_t i = 0; i  < kept_kp_count; ++i)
      pt_to_track[i] = features::PointFeature(m_nextKeypoints[i].pt.x, m_nextKeypoints[i].pt.y);

    return kept_kp_count != 0;
    // Return false if no point can be added
  }
};

} // namespace VO
} // namespace openMVG


#endif // TRACKER_OPENCV_KLT_HPP
