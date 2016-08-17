#pragma once

#include <openMVG/dataio/FeedProvider.hpp>

#include <opencv2/opencv.hpp>

#include <vector>
#include <string>

namespace openMVG{
namespace calibration{

/**
 * @brief This function exports undistorted images.
 *
 * @param[in] feed
 * @param[in] debugFolder
 * @param[in] exportFrames
 * @param[in] cameraMatrix
 * @param[in] distCoeffs
 * @param[in] imageSize
 * @param[in] suffix
 */
void exportImages(openMVG::dataio::FeedProvider& feed,
                  const std::string& debugFolder,
                  const std::vector<std::size_t>& exportFrames,
                  const cv::Mat& cameraMatrix,
                  const cv::Mat& distCoeffs,
                  const cv::Size& imageSize,
                  const std::string& suffix = "_undistort.png");

/**
 * @brief This debug function lets the user to see the undistorted images.
 *
 * @param[in] debugSelectedImgFolder
 * @param[in] debugRejectedImgFolder
 * @param[in] feed
 * @param[in] calibInputFrames
 * @param[in] rejectInputFrames
 * @param[in] remainingImagesIndexes
 * @param[in] cameraMatrix
 * @param[in] distCoeffs
 * @param[in] imageSize
 */
void exportDebug(const std::string& debugSelectedImgFolder,
                 const std::string& debugRejectedImgFolder,
                 openMVG::dataio::FeedProvider& feed,
                 const std::vector<std::size_t>& calibInputFrames,
                 const std::vector<std::size_t>& rejectInputFrames,
                 const std::vector<std::size_t>& remainingImagesIndexes,
                 const cv::Mat& cameraMatrix,
                 const cv::Mat& distCoeffs,
                 const cv::Size& imageSize);

/**
 * @brief This function saves the parameters' camera into a txt file.
 *
 * @param[out] filename
 * @param[in] imageSize
 * @param[in] cameraMatrix
 * @param[in] distCoeffs
 */
void saveCameraParamsToPlainTxt(const std::string& filename,
                                const cv::Size& imageSize,
                                const cv::Mat& cameraMatrix,
                                const cv::Mat& distCoeffs);

/**
 * @brief This function saves some parameters' camera into a txt file.
 *
 * @param[in] filename
 * @param[in] imageSize
 * @param[in] boardSize
 * @param[in] squareSize
 * @param[in] aspectRatio
 * @param[in] cvCalibFlags
 * @param[in] cameraMatrix
 * @param[in] distCoeffs
 * @param[in] rvecs
 * @param[in] tvecs
 * @param[in] reprojErrs
 * @param[in] imagePoints
 * @param[in] totalAvgErr
 */
void saveCameraParams(const std::string& filename,
                      const cv::Size& imageSize, const cv::Size& boardSize,
                      float squareSize, float aspectRatio, int cvCalibFlags,
                      const cv::Mat& cameraMatrix, const cv::Mat& distCoeffs,
                      const std::vector<cv::Mat>& rvecs, const std::vector<cv::Mat>& tvecs,
                      const std::vector<float>& reprojErrs,
                      const std::vector<std::vector<cv::Point2f> >& imagePoints,
                      double totalAvgErr);

}//namespace calibration
}//namespace openMVG

