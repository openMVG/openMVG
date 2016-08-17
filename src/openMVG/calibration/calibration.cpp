#include "openMVG/calibration/calibration.hpp"

#include <ctime>
#include <algorithm>
#include <iostream>

namespace openMVG{
namespace calibration{

double computeReprojectionErrors(const std::vector<std::vector<cv::Point3f> >& objectPoints,
                                 const std::vector<std::vector<cv::Point2f> >& imagePoints,
                                 const std::vector<cv::Mat>& rvecs, const std::vector<cv::Mat>& tvecs,
                                 const cv::Mat& cameraMatrix, const cv::Mat& distCoeffs,
                                 std::vector<float>& perViewErrors)
{
  std::vector<cv::Point2f> imagePoints2;
  int totalPoints = 0;
  double totalErr = 0;
  perViewErrors.resize(objectPoints.size());

  for (std::size_t i = 0; i < (int) objectPoints.size(); i++)
  {
    cv::projectPoints(cv::Mat(objectPoints[i]), rvecs[i], tvecs[i],
                      cameraMatrix, distCoeffs, imagePoints2);
    const double err = cv::norm(cv::Mat(imagePoints[i]), cv::Mat(imagePoints2), CV_L2);
    const std::size_t n = (int) objectPoints[i].size();
    perViewErrors[i] = (float) std::sqrt(err * err / n);
    totalErr += err*err;
    totalPoints += n;
  }
  return std::sqrt(totalErr / totalPoints);
}

bool runCalibration(const std::vector<std::vector<cv::Point2f> >& imagePoints,
                    const std::vector<std::vector<cv::Point3f> >& objectPoints,
                    const cv::Size imageSize,
                    float aspectRatio,
                    int cvCalibFlags,
                    cv::Mat& cameraMatrix,
                    cv::Mat& distCoeffs,
                    std::vector<cv::Mat>& rvecs,
                    std::vector<cv::Mat>& tvecs,
                    std::vector<float>& reprojErrs,
                    double& totalAvgErr)
{
  rvecs.resize(0);
  tvecs.resize(0);
  reprojErrs.resize(0);
  cameraMatrix = cv::Mat::eye(3, 3, CV_64F);
  if (cvCalibFlags & CV_CALIB_FIX_ASPECT_RATIO)
    cameraMatrix.at<double>(0, 0) = aspectRatio;

  distCoeffs = cv::Mat::zeros(8, 1, CV_64F);

  std::clock_t startrC = std::clock();

  const double rms = cv::calibrateCamera(objectPoints, imagePoints, imageSize, cameraMatrix,
                                         distCoeffs, rvecs, tvecs, cvCalibFlags | CV_CALIB_FIX_K4 | CV_CALIB_FIX_K5 | CV_CALIB_FIX_K6);
  std::clock_t durationrC = (std::clock() - startrC) / (double) CLOCKS_PER_SEC;
  std::cout << "  calibrateCamera duration: " << durationrC << std::endl;

  printf("RMS error reported by calibrateCamera: %g\n", rms);
  bool ok = cv::checkRange(cameraMatrix) && cv::checkRange(distCoeffs);

  startrC = std::clock();

  totalAvgErr = computeReprojectionErrors(objectPoints, imagePoints,
                                          rvecs, tvecs, cameraMatrix, distCoeffs, reprojErrs);
  durationrC = (std::clock() - startrC) / (double) CLOCKS_PER_SEC;
  std::cout << "  computeReprojectionErrors duration: " << durationrC << std::endl;

  return ok;
}

bool calibrationIterativeOptimization(std::vector<std::vector<cv::Point2f> >& calibImagePoints,
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
                        std::vector<std::size_t>& rejectInputFrames)
{
  std::size_t calibIteration = 0;
  bool calibSucceeded = false;
  do
  {
    // Estimate the camera calibration
    std::cout << "Calibration iteration " << calibIteration << " with " << calibImagePoints.size() << " frames." << std::endl;
    calibSucceeded = runCalibration(calibImagePoints, calibObjectPoints, imageSize,
                                    aspectRatio, cvCalibFlags, cameraMatrix, distCoeffs,
                                    rvecs, tvecs, reprojErrs, totalAvgErr);

    if (totalAvgErr <= maxTotalAvgErr)
    {
      std::cout << "The calibration succeed with an average error that respects the maxTotalAvgErr." << std::endl;
      break;
    }
    else if (calibInputFrames.size() < minInputFrames)
    {
      std::cout << "Not enough valid input image (" << calibInputFrames.size() << ") to continue the refinement." << std::endl;
      break;
    }
    else if (calibSucceeded)
    {
      // Filter the successfully calibrated images to keep the best ones
      // in order to refine the calibration.
      // For instance, remove blurry images which introduce imprecision.

      std::vector<float> globalScores;
      for (std::size_t i = 0; i < calibInputFrames.size(); ++i)
      {
        globalScores.push_back(reprojErrs[i] * calibImageScore[i]);
      }

      const auto minMaxError = std::minmax_element(globalScores.begin(), globalScores.end());
      std::cout << "minMaxError: " << *minMaxError.first << ", " << *minMaxError.second << std::endl;
      if (*minMaxError.first == *minMaxError.second)
      {
        std::cout << "Same error on all images: " << *minMaxError.first << std::endl;
        for (float f : globalScores)
          std::cout << "f: " << f << std::endl;
        break;
      }
      // We only keep the frames with N% of the largest error.
      const float errorThreshold = *minMaxError.first + 0.8 * (*minMaxError.second - *minMaxError.first);
      std::vector<std::vector<cv::Point2f> > filteredImagePoints;
      std::vector<std::vector<cv::Point3f> > filteredObjectPoints;
      std::vector<std::size_t> filteredInputFrames;
      std::vector<std::size_t> tmpRejectInputFrames;
      std::vector<float> filteredImageScores;

      for (std::size_t i = 0; i < calibImagePoints.size(); ++i)
      {
        if (globalScores[i] < errorThreshold)
        {
          filteredImagePoints.push_back(calibImagePoints[i]);
          filteredObjectPoints.push_back(calibObjectPoints[i]);
          filteredInputFrames.push_back(calibInputFrames[i]);
          filteredImageScores.push_back(calibImageScore[i]);
        }
        else
        {
          // We collect rejected frames for debug purpose
          tmpRejectInputFrames.push_back(calibInputFrames[i]);
        }
      }
      if (filteredImagePoints.size() < minInputFrames)
      {
        std::cout << "Not enough filtered input images (filtered: " << filteredImagePoints.size() << ", rejected:" << tmpRejectInputFrames.size() << ") to continue the refinement." << std::endl;
        break;
      }
      if (calibImagePoints.size() == filteredImagePoints.size())
      {
        // Convergence reached
        std::cout << "Convergence reached." << std::endl;
        break;
      }
      calibImagePoints.swap(filteredImagePoints);
      calibObjectPoints.swap(filteredObjectPoints);
      calibInputFrames.swap(filteredInputFrames);
      calibImageScore.swap(filteredImageScores);
      rejectInputFrames.insert(rejectInputFrames.end(), tmpRejectInputFrames.begin(), tmpRejectInputFrames.end());
    }
    ++calibIteration;
  }
  while (calibSucceeded);

  std::cout << "Calibration done with " << calibIteration << " iterations." << std::endl;
  std::cout << "Average reprojection error is " << totalAvgErr << std::endl;
  std::cout << (calibSucceeded ? "Calibration succeeded" : "Calibration failed") << std::endl;

  return calibSucceeded;
}

}//namespace calibration
}//namespace openMVG
