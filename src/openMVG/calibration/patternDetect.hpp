#pragma once

#include <vector>
#include <iostream>
#include <opencv2/opencv.hpp>

namespace openMVG{
namespace calibration{

enum Pattern
{
  CHESSBOARD = 0,
  CIRCLES_GRID,
  ASYMMETRIC_CIRCLES_GRID
#ifdef HAVE_CCTAG
  , CCTAG_GRID
#endif
};

/**
 * @brief Read pattern from console.
 * 
 * @param[in,out] stream
 * @param[in] pattern
 * @return stream
 */
std::istream& operator>>(std::istream &stream, Pattern &pattern);

/**
 * @brief Write pattern to console.
 * 
 * @param[out] stream
 * @param[in] pattern
 * @return stream
 */
std::ostream& operator<<(std::ostream &stream, const Pattern pattern);

/**
 * @brief This function detects the checkerboard in images
 *
 * @param[in] pattern
 * @param[in] viewGray
 * @param[in] boardSize
 * @param[out] pointbuf
 * @return True if the pattern is found, otherwise false.
 */
bool findPattern(const Pattern& pattern, const cv::Mat& viewGray, const cv::Size& boardSize, std::vector<cv::Point2f>& pointbuf);

/**
 * @brief This function computes the points' coordinates of the checkerboard.
 *
 * @param[in] boardSize
 * @param[in] squareSize
 * @param[out] corners
 * @param[in] pattern
 */
void calcChessboardCorners(const cv::Size& boardSize, const float& squareSize,
                           std::vector<cv::Point3f>& corners, Pattern pattern);

/**
 * @brief This function creates an object which stores all the points of the images.
 *
 * @param[in] boardSize
 * @param[in] pattern
 * @param[in] squareSize
 * @param[in] imagePoints
 * @param[out] objectPoints
 */
void computeObjectPoints(const cv::Size& boardSize,
                         Pattern pattern,
                         const float& squareSize,
                         const std::vector<std::vector<cv::Point2f> >& imagePoints,
                         std::vector<std::vector<cv::Point3f> >& objectPoints);

}//namespace calibration
}//namespace openMVG


