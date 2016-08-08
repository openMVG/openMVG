#pragma once

#include <vector>
#include <opencv2/opencv.hpp>

namespace openMVG {
namespace calibration {

/**
 * @brief This function computes the average of the reprojection errors.
 *
 * @param[in] objectPoints
 * @param[in] imagePoints
 * @param[in] rvecs
 * @param[in] tvecs
 * @param[in] cameraMatrix
 * @param[in] distCoeffs
 * @param[out] perViewErrors
 * @return The average of the reprojection errors.
 */
static double computeReprojectionErrors(
                                        const std::vector<std::vector<cv::Point3f> >& objectPoints,
                                        const std::vector<std::vector<cv::Point2f> >& imagePoints,
                                        const std::vector<cv::Mat>& rvecs, const std::vector<cv::Mat>& tvecs,
                                        const cv::Mat& cameraMatrix, const cv::Mat& distCoeffs,
                                        std::vector<float>& perViewErrors);

/**
 * @brief This function calibrates the camera.
 *
 * @param[in] imagePoints
 * @param[in] objectPoints
 * @param[in] imageSize
 * @param[in] aspectRatio
 * @param[in] cvCalibFlags
 * @param[in] cameraMatrix
 * @param[in] distCoeffs
 * @param[in] rvecs
 * @param[in] tvecs
 * @param[in] reprojErrs
 * @param[out] totalAvgErr
 * @return True if the calibration is a success, otherwise false.
 */
static bool runCalibration(const std::vector<std::vector<cv::Point2f> >& imagePoints,
                           const std::vector<std::vector<cv::Point3f> >& objectPoints,
                           cv::Size imageSize,
                           float aspectRatio,
                           int cvCalibFlags, cv::Mat& cameraMatrix,
                           cv::Mat& distCoeffs,
                           std::vector<cv::Mat>& rvecs,
                           std::vector<cv::Mat>& tvecs,
                           std::vector<float>& reprojErrs,
                           double& totalAvgErr);

/**
 * @brief This function is the refinement loop of the calibration.
 *
 * @param[in,out] calibImagePoints
 * @param[in,out] calibObjectPoints
 * @param[in] imageSize
 * @param[in] aspectRatio
 * @param[in] cvCalibFlags
 * @param[in] cameraMatrix
 * @param[in] distCoeffs
 * @param[in] rvecs
 * @param[in] tvecs
 * @param[in] reprojErrs
 * @param[in] totalAvgErr
 * @param[in] maxTotalAvgErr
 * @param[in] minInputFrames
 * @param[in,out] calibInputFrames
 * @param[in,out] calibImageScore
 * @param[out] rejectInputFrames
 * @return True if the calibration is a success, otherwise false.
 */
int calibrationIterativeOptimization(std::vector<std::vector<cv::Point2f> >& calibImagePoints,
                        std::vector<std::vector<cv::Point3f> >& calibObjectPoints,
                        const cv::Size& imageSize,
                        float aspectRatio,
                        int cvCalibFlags,
                        cv::Mat& cameraMatrix,
                        cv::Mat& distCoeffs,
                        std::vector<cv::Mat>& rvecs,
                        std::vector<cv::Mat>& tvecs,
                        std::vector<float>& reprojErrs,
                        double& totalAvgErr,
                        const double& maxTotalAvgErr,
                        const std::size_t& minInputFrames,
                        std::vector<std::size_t>& calibInputFrames,
                        std::vector<float>& calibImageScore,
                        std::vector<std::size_t>& rejectInputFrames);

}//namespace calibration
}//namespace openMVG
