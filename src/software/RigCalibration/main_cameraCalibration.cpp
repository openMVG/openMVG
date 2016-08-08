#include <openMVG/dataio/FeedProvider.hpp>
#include <openMVG/cameras/Camera_undistort_image.hpp>
#include <openMVG/cameras/Camera_Pinhole_Radial.hpp>
#include <openMVG/image/image_io.hpp>
#include <openMVG/calibration/patternDetect.hpp>
#include <openMVG/calibration/bestImages.hpp>
#include <openMVG/calibration/calibration.hpp>
#include <openMVG/calibration/exportData.hpp>

#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/algorithm/string/case_conv.hpp>

#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/core/eigen.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/calib3d/calib3d.hpp>

#include <stdio.h>
#include <ctime>
#include <cstdio>
#include <string>
#include <cctype>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <exception>
#include <map>
#include <limits>

namespace bfs = boost::filesystem;
namespace po = boost::program_options;

///**
// * @brief This function exports undistorted images.
// *
// * @param[in] feed
// * @param[in] debugFolder
// * @param[in] exportFrames
// * @param[in] cameraMatrix
// * @param[in] distCoeffs
// * @param[in] imageSize
// * @param[in] suffix
// */
//void exportImages(openMVG::dataio::FeedProvider& feed,
//                  const std::string& debugFolder,
//                  const std::vector<std::size_t>& exportFrames,
//                  const cv::Mat& cameraMatrix,
//                  const cv::Mat& distCoeffs,
//                  const cv::Size& imageSize,
//                  const std::string suffix = "_undistort.png")
//{
//  std::vector<int> export_params;
//  openMVG::image::Image<unsigned char> inputImage;
//  openMVG::image::Image<unsigned char> outputImage;
//  std::string currentImgName;
//  openMVG::cameras::Pinhole_Intrinsic_Radial_K3 queryIntrinsics;
//  bool hasIntrinsics;
//
//  export_params.push_back(CV_IMWRITE_JPEG_QUALITY);
//  export_params.push_back(100);
//
//  openMVG::cameras::Pinhole_Intrinsic_Radial_K3 camera(
//                                                       imageSize.width, imageSize.height,
//                                                       cameraMatrix.at<double>(0, 0), cameraMatrix.at<double>(0, 2), cameraMatrix.at<double>(1, 2),
//                                                       distCoeffs.at<double>(0), distCoeffs.at<double>(1), distCoeffs.at<double>(4));
//  std::cout << "Coefficients matrix :\n " << distCoeffs << std::endl;
//  std::cout << "Exporting images ..." << std::endl;
//  for (std::size_t currentFrame : exportFrames)
//  {
//    feed.goToFrame(currentFrame);
//    feed.readImage(inputImage, queryIntrinsics, currentImgName, hasIntrinsics);
//
//    // drawChessboardCorners(view, boardSize, cv::Mat(pointbuf), found);
//
//    openMVG::cameras::UndistortImage(inputImage, &camera, outputImage, openMVG::image::BLACK);
//    const bfs::path imagePath = bfs::path(debugFolder) / (std::to_string(currentFrame) + suffix);
//    const bool exportStatus = openMVG::image::WriteImage(imagePath.c_str(), outputImage);
//    if (!exportStatus)
//    {
//      std::cerr << "Failed to export: " << imagePath << std::endl;
//    }
//  }
//  std::cout << "... finished" << std::endl;
//}
//
///**
// * @brief This debug function lets the user to see the undistorted images.
// *
// * @param[in] debugSelectedImgFolder
// * @param[in] debugRejectedImgFolder
// * @param[in] feed
// * @param[in] calibInputFrames
// * @param[in] rejectInputFrames
// * @param[in] remainingImagesIndexes
// * @param[in] cameraMatrix
// * @param[in] distCoeffs
// * @param[in] imageSize
// */
//void exportDebug(const std::string& debugSelectedImgFolder,
//                 const std::string& debugRejectedImgFolder,
//                 openMVG::dataio::FeedProvider& feed,
//                 const std::vector<std::size_t>& calibInputFrames,
//                 const std::vector<std::size_t>& rejectInputFrames,
//                 const std::vector<std::size_t>& remainingImagesIndexes,
//                 const cv::Mat& cameraMatrix,
//                 const cv::Mat& distCoeffs,
//                 const cv::Size& imageSize)
//{
//  std::clock_t startDebug = std::clock();
//  double durationDebug;
//
//  if (!debugSelectedImgFolder.empty())
//  {
//
//    startDebug = std::clock();
//    exportImages(feed, debugSelectedImgFolder, calibInputFrames,
//                 cameraMatrix, distCoeffs, imageSize, "_undistort.png");
//    durationDebug = (std::clock() - startDebug) / (double) CLOCKS_PER_SEC;
//    std::cout << "Export debug of selected frames, duration: " << durationDebug << std::endl;
//  }
//
//  if (!debugRejectedImgFolder.empty())
//  {
//    startDebug = std::clock();
//    exportImages(feed, debugRejectedImgFolder, rejectInputFrames,
//                 cameraMatrix, distCoeffs, imageSize, "_rejected_undistort.png");
//    durationDebug = (std::clock() - startDebug) / (double) CLOCKS_PER_SEC;
//    std::cout << "Export debug of rejected frames, duration: " << durationDebug << std::endl;
//  }
//
//  if (!debugRejectedImgFolder.empty())
//  {
//    startDebug = std::clock();
//    exportImages(feed, debugRejectedImgFolder, remainingImagesIndexes,
//                 cameraMatrix, distCoeffs, imageSize, "_not_selected_undistort.png");
//    durationDebug = (std::clock() - startDebug) / (double) CLOCKS_PER_SEC;
//    std::cout << "Export debug of not selected frames, duration: " << durationDebug << std::endl;
//  }
//}

///**
// * @brief This function computes the average of the reprojection errors.
// *
// * @param[in] objectPoints
// * @param[in] imagePoints
// * @param[in] rvecs
// * @param[in] tvecs
// * @param[in] cameraMatrix
// * @param[in] distCoeffs
// * @param[out] perViewErrors
// * @return The average of the reprojection errors.
// */
//static double computeReprojectionErrors(
//                                        const std::vector<std::vector<cv::Point3f> >& objectPoints,
//                                        const std::vector<std::vector<cv::Point2f> >& imagePoints,
//                                        const std::vector<cv::Mat>& rvecs, const std::vector<cv::Mat>& tvecs,
//                                        const cv::Mat& cameraMatrix, const cv::Mat& distCoeffs,
//                                        std::vector<float>& perViewErrors)
//{
//  std::vector<cv::Point2f> imagePoints2;
//  int i, totalPoints = 0;
//  double err, totalErr = 0;
//  perViewErrors.resize(objectPoints.size());
//
//  for (i = 0; i < (int) objectPoints.size(); i++)
//  {
//    cv::projectPoints(cv::Mat(objectPoints[i]), rvecs[i], tvecs[i],
//                      cameraMatrix, distCoeffs, imagePoints2);
//    err = cv::norm(cv::Mat(imagePoints[i]), cv::Mat(imagePoints2), CV_L2);
//    int n = (int) objectPoints[i].size();
//    perViewErrors[i] = (float) std::sqrt(err * err / n);
//    totalErr += err*err;
//    totalPoints += n;
//  }
//
//  return std::sqrt(totalErr / totalPoints);
//}

///**
// * @brief This function computes the points' coordinates of the checkerboard.
// *
// * @param[in] boardSize
// * @param[in] squareSize
// * @param[out] corners
// * @param[in] pattern
// */
//static void calcChessboardCorners(cv::Size boardSize, float squareSize,
//                                  std::vector<cv::Point3f>& corners, Pattern pattern = CHESSBOARD)
//{
//  corners.resize(0);
//
//  switch (pattern)
//  {
//  case CHESSBOARD:
//  case CIRCLES_GRID:
//    for (int i = 0; i < boardSize.height; i++)
//      for (int j = 0; j < boardSize.width; j++)
//        corners.push_back(cv::Point3f(float(j * squareSize),
//                                      float(i * squareSize), 0));
//    break;
//
//  case ASYMMETRIC_CIRCLES_GRID:
//    for (int i = 0; i < boardSize.height; i++)
//      for (int j = 0; j < boardSize.width; j++)
//        corners.push_back(cv::Point3f(float((2 * j + i % 2) * squareSize),
//                                      float(i * squareSize), 0));
//    break;
//
//#ifdef HAVE_CCTAG
//  case CCTAG_GRID:
//    throw std::invalid_argument("CCTag not implemented.");
//    break;
//#endif
//
//  default:
//    throw std::invalid_argument("Unknown pattern type.");
//  }
//}

///**
// * @brief This function detects the checkerboard in images
// *
// * @param[in] pattern
// * @param[out] found
// * @param[in] viewGray
// * @param[in] boardSize
// * @param[in] pointbuf
// * @return True if the pattern is found, otherwise false.
// */
//int findPattern(const Pattern& pattern, bool& found, const cv::Mat& viewGray, const cv::Size& boardSize, std::vector<cv::Point2f>& pointbuf)
//{
//  std::clock_t startCh;
//  double durationCh;
//
//  switch (pattern)
//  {
//  case CHESSBOARD:
//    startCh = std::clock();
//
//    found = cv::findChessboardCorners(viewGray, boardSize, pointbuf,
//                                      CV_CALIB_CB_ADAPTIVE_THRESH | CV_CALIB_CB_FAST_CHECK | CV_CALIB_CB_NORMALIZE_IMAGE);
//    durationCh = (std::clock() - startCh) / (double) CLOCKS_PER_SEC;
//    std::cout << "Find chessboard corners' duration: " << durationCh << std::endl;
//    startCh = std::clock();
//
//    // improve the found corners' coordinate accuracy
//    if (found)
//      cv::cornerSubPix(viewGray, pointbuf, cv::Size(11, 11), cv::Size(-1, -1),
//                       cv::TermCriteria(CV_TERMCRIT_EPS + CV_TERMCRIT_ITER, 30, 0.1));
//
//    durationCh = (std::clock() - startCh) / (double) CLOCKS_PER_SEC;
//    std::cout << "Refine chessboard corners' duration: " << durationCh << std::endl;
//    break;
//
//  case CIRCLES_GRID:
//    startCh = std::clock();
//
//    found = cv::findCirclesGrid(viewGray, boardSize, pointbuf);
//
//    durationCh = (std::clock() - startCh) / (double) CLOCKS_PER_SEC;
//    std::cout << "Find circles grid duration: " << durationCh << std::endl;
//    break;
//
//  case ASYMMETRIC_CIRCLES_GRID:
//    startCh = std::clock();
//
//    found = cv::findCirclesGrid(viewGray, boardSize, pointbuf, cv::CALIB_CB_ASYMMETRIC_GRID);
//
//    durationCh = (std::clock() - startCh) / (double) CLOCKS_PER_SEC;
//    std::cout << "Find asymmetric circles grid duration: " << durationCh << std::endl;
//    break;
//
//#ifdef HAVE_CCTAG
//  case CCTAG_GRID:
//    throw std::invalid_argument("CCTag calibration not implemented.");
//    break;
//#endif
//
//  default:
//    std::cerr << "Unknown pattern type" << std::endl;
//    return EXIT_FAILURE;
//  }
//  return EXIT_SUCCESS;
//}

///**
// * @brief This function creates an object which stores all the points of the images.
// *
// * @param[in] boardSize
// * @param[in] pattern
// * @param[in] squareSize
// * @param[in] imagePoints
// * @param[out] objectPoints
// */
//static void computeObjectPoints(
//                                cv::Size boardSize,
//                                Pattern pattern,
//                                float squareSize,
//                                const std::vector<std::vector<cv::Point2f> >& imagePoints,
//                                std::vector<std::vector<cv::Point3f> >& objectPoints)
//{
//  std::vector<cv::Point3f> templateObjectPoints;
//
//  // Generate the object points coordinates
//  calcChessboardCorners(boardSize, squareSize, templateObjectPoints, pattern);
//
//  // Assign the template to all items
//  objectPoints.resize(imagePoints.size(), templateObjectPoints);
//}

///**
// * @brief This function calibrates the camera.
// *
// * @param[in] imagePoints
// * @param[in] objectPoints
// * @param[in] imageSize
// * @param[in] aspectRatio
// * @param[in] cvCalibFlags
// * @param[in] cameraMatrix
// * @param[in] distCoeffs
// * @param[in] rvecs
// * @param[in] tvecs
// * @param[in] reprojErrs
// * @param[out] totalAvgErr
// * @return True if the calibration is a success, otherwise false.
// */
//static bool runCalibration(const std::vector<std::vector<cv::Point2f> >& imagePoints,
//                           const std::vector<std::vector<cv::Point3f> >& objectPoints,
//                           cv::Size imageSize,
//                           float aspectRatio,
//                           int cvCalibFlags, cv::Mat& cameraMatrix,
//                           cv::Mat& distCoeffs,
//                           std::vector<cv::Mat>& rvecs,
//                           std::vector<cv::Mat>& tvecs,
//                           std::vector<float>& reprojErrs,
//                           double& totalAvgErr)
//{
//  rvecs.resize(0);
//  tvecs.resize(0);
//  reprojErrs.resize(0);
//  cameraMatrix = cv::Mat::eye(3, 3, CV_64F);
//  if (cvCalibFlags & CV_CALIB_FIX_ASPECT_RATIO)
//    cameraMatrix.at<double>(0, 0) = aspectRatio;
//
//  distCoeffs = cv::Mat::zeros(8, 1, CV_64F);
//
//  std::clock_t startrC = std::clock();
//
//  const double rms = cv::calibrateCamera(objectPoints, imagePoints, imageSize, cameraMatrix,
//                                         distCoeffs, rvecs, tvecs, cvCalibFlags | CV_CALIB_FIX_K4 | CV_CALIB_FIX_K5 | CV_CALIB_FIX_K6);
//  std::clock_t durationrC = (std::clock() - startrC) / (double) CLOCKS_PER_SEC;
//  std::cout << "  calibrateCamera duration: " << durationrC << std::endl;
//
//  printf("RMS error reported by calibrateCamera: %g\n", rms);
//  bool ok = cv::checkRange(cameraMatrix) && cv::checkRange(distCoeffs);
//
//  startrC = std::clock();
//
//  totalAvgErr = computeReprojectionErrors(objectPoints, imagePoints,
//                                          rvecs, tvecs, cameraMatrix, distCoeffs, reprojErrs);
//  durationrC = (std::clock() - startrC) / (double) CLOCKS_PER_SEC;
//  std::cout << "  computeReprojectionErrors duration: " << durationrC << std::endl;
//
//  return ok;
//}

///**
// * @brief This function is the refinement loop of the calibration.
// *
// * @param[in,out] calibImagePoints
// * @param[in,out] calibObjectPoints
// * @param[in] imageSize
// * @param[in] aspectRatio
// * @param[in] cvCalibFlags
// * @param[in] cameraMatrix
// * @param[in] distCoeffs
// * @param[in] rvecs
// * @param[in] tvecs
// * @param[in] reprojErrs
// * @param[in] totalAvgErr
// * @param[in] maxTotalAvgErr
// * @param[in] minInputFrames
// * @param[in,out] calibInputFrames
// * @param[in,out] calibImageScore
// * @param[out] rejectInputFrames
// * @return True if the calibration is a success, otherwise false.
// */
//int calibrationIterativeOptimization(std::vector<std::vector<cv::Point2f> >& calibImagePoints,
//                        std::vector<std::vector<cv::Point3f> >& calibObjectPoints,
//                        const cv::Size& imageSize,
//                        float aspectRatio,
//                        int cvCalibFlags,
//                        cv::Mat& cameraMatrix,
//                        cv::Mat& distCoeffs,
//                        std::vector<cv::Mat>& rvecs,
//                        std::vector<cv::Mat>& tvecs,
//                        std::vector<float>& reprojErrs,
//                        double& totalAvgErr,
//                        const double& maxTotalAvgErr,
//                        const std::size_t& minInputFrames,
//                        std::vector<std::size_t>& calibInputFrames,
//                        std::vector<float>& calibImageScore,
//                        std::vector<std::size_t>& rejectInputFrames)
//{
//  std::size_t calibIteration = 0;
//  bool calibSucceeded = false;
//  do
//  {
//    // Estimate the camera calibration
//    std::cout << "Calibration iteration " << calibIteration << " with " << calibImagePoints.size() << " frames." << std::endl;
//    calibSucceeded = runCalibration(calibImagePoints, calibObjectPoints, imageSize,
//                                    aspectRatio, cvCalibFlags, cameraMatrix, distCoeffs,
//                                    rvecs, tvecs, reprojErrs, totalAvgErr);
//
//    if (totalAvgErr <= maxTotalAvgErr)
//    {
//      std::cout << "The calibration succeed with an average error that respects the maxTotalAvgErr." << std::endl;
//      break;
//    }
//    else if (calibInputFrames.size() < minInputFrames)
//    {
//      std::cout << "Not enough valid input image (" << calibInputFrames.size() << ") to continue the refinement." << std::endl;
//      break;
//    }
//    else if (calibSucceeded)
//    {
//      // Filter the successfully calibrated images to keep the best ones
//      // in order to refine the calibration.
//      // For instance, remove blurry images which introduce imprecision.
//
//      std::vector<float> globalScores;
//      for (int i = 0; i < calibInputFrames.size(); ++i)
//      {
//        globalScores.push_back(reprojErrs[i] * calibImageScore[i]);
//      }
//
//      const auto minMaxError = std::minmax_element(globalScores.begin(), globalScores.end());
//      std::cout << "minMaxError: " << *minMaxError.first << ", " << *minMaxError.second << std::endl;
//      if (*minMaxError.first == *minMaxError.second)
//      {
//        std::cout << "Same error on all images: " << *minMaxError.first << std::endl;
//        for (float f : globalScores)
//          std::cout << "f: " << f << std::endl;
//        break;
//      }
//      // We only keep the frames with N% of the largest error.
//      const float errorThreshold = *minMaxError.first + 0.8 * (*minMaxError.second - *minMaxError.first);
//      std::vector<std::vector<cv::Point2f> > filteredImagePoints;
//      std::vector<std::vector<cv::Point3f> > filteredObjectPoints;
//      std::vector<std::size_t> filteredInputFrames;
//      std::vector<std::size_t> tmpRejectInputFrames;
//      std::vector<float> filteredImageScores;
//
//      for (std::size_t i = 0; i < calibImagePoints.size(); ++i)
//      {
//        if (globalScores[i] < errorThreshold)
//        {
//          filteredImagePoints.push_back(calibImagePoints[i]);
//          filteredObjectPoints.push_back(calibObjectPoints[i]);
//          filteredInputFrames.push_back(calibInputFrames[i]);
//          filteredImageScores.push_back(calibImageScore[i]);
//        }
//        else
//        {
//          // We collect rejected frames for debug purpose
//          tmpRejectInputFrames.push_back(calibInputFrames[i]);
//        }
//      }
//      if (filteredImagePoints.size() < minInputFrames)
//      {
//        std::cout << "Not enough filtered input images (filtered: " << filteredImagePoints.size() << ", rejected:" << tmpRejectInputFrames.size() << ") to continue the refinement." << std::endl;
//        break;
//      }
//      if (calibImagePoints.size() == filteredImagePoints.size())
//      {
//        // Convergence reached
//        std::cout << "Convergence reached." << std::endl;
//        break;
//      }
//      calibImagePoints.swap(filteredImagePoints);
//      calibObjectPoints.swap(filteredObjectPoints);
//      calibInputFrames.swap(filteredInputFrames);
//      calibImageScore.swap(filteredImageScores);
//      rejectInputFrames.insert(rejectInputFrames.end(), tmpRejectInputFrames.begin(), tmpRejectInputFrames.end());
//    }
//    ++calibIteration;
//  }
//  while (calibSucceeded);
//
//  std::cout << "Calibration done with " << calibIteration << " iterations." << std::endl;
//  std::cout << "Average reprojection error is " << totalAvgErr << std::endl;
//  std::cout << (calibSucceeded ? "Calibration succeeded" : "Calibration failed") << std::endl;
//
//  if (!calibSucceeded)
//    return EXIT_FAILURE;
//
//  return EXIT_SUCCESS;
//}

///**
// * @brief This function saves the parameters' camera into a txt file.
// *
// * @param[out] filename
// * @param[in] imageSize
// * @param[in] cameraMatrix
// * @param[in] distCoeffs
// */
//static void saveCameraParamsToPlainTxt(const std::string &filename,
//                                       const cv::Size &imageSize,
//                                       const cv::Mat& cameraMatrix,
//                                       const cv::Mat& distCoeffs)
//{
//  std::ofstream fs(filename, std::ios::out);
//  if (!fs.is_open())
//  {
//    std::cerr << "Unable to create the calibration file " << filename << std::endl;
//    throw std::invalid_argument("Unable to create the calibration file " + filename);
//  }
//
//  // the structure of the file is
//  // int #image width
//  // int #image height
//  // double #focal
//  // double #ppx principal point x-coord
//  // double #ppy principal point y-coord
//  // double #k0
//  // double #k1
//  // double #k2
//  fs << imageSize.width << std::endl;
//  fs << imageSize.height << std::endl;
//  if (cameraMatrix.type() == cv::DataType<double>::type)
//  {
//    fs << (cameraMatrix.at<double>(0, 0) + cameraMatrix.at<double>(1, 1)) / 2 << std::endl;
//    fs << cameraMatrix.at<double>(0, 2) << std::endl;
//    fs << cameraMatrix.at<double>(1, 2) << std::endl;
//  }
//  else
//  {
//    fs << (cameraMatrix.at<float>(0, 0) + cameraMatrix.at<float>(1, 1)) / 2 << std::endl;
//    fs << cameraMatrix.at<float>(0, 2) << std::endl;
//    fs << cameraMatrix.at<float>(1, 2) << std::endl;
//  }
//  if (distCoeffs.type() == cv::DataType<double>::type)
//  {
//    fs << distCoeffs.at<double>(0) << std::endl;
//    fs << distCoeffs.at<double>(1) << std::endl;
//    fs << distCoeffs.at<double>(4) << std::endl;
//  }
//  else
//  {
//    fs << distCoeffs.at<float>(0) << std::endl;
//    fs << distCoeffs.at<float>(1) << std::endl;
//    fs << distCoeffs.at<float>(4) << std::endl;
//  }
//  fs.close();
//}
///**
// * @brief This function saves some parameters' camera into a txt file.
// *
// * @param[in] filename
// * @param[in] imageSize
// * @param[in] boardSize
// * @param[in] squareSize
// * @param[in] aspectRatio
// * @param[in] cvCalibFlags
// * @param[in] cameraMatrix
// * @param[in] distCoeffs
// * @param[in] rvecs
// * @param[in] tvecs
// * @param[in] reprojErrs
// * @param[in] imagePoints
// * @param[in] totalAvgErr
// */
//static void saveCameraParams(const std::string& filename,
//                             cv::Size imageSize, cv::Size boardSize,
//                             float squareSize, float aspectRatio, int cvCalibFlags,
//                             const cv::Mat& cameraMatrix, const cv::Mat& distCoeffs,
//                             const std::vector<cv::Mat>& rvecs, const std::vector<cv::Mat>& tvecs,
//                             const std::vector<float>& reprojErrs,
//                             const std::vector<std::vector<cv::Point2f> >& imagePoints,
//                             double totalAvgErr)
//{
//  cv::FileStorage fs(filename, cv::FileStorage::WRITE);
//
//  time_t tt;
//  time(&tt);
//  struct tm *t2 = localtime(&tt);
//  char buf[1024];
//  strftime(buf, sizeof (buf) - 1, "%c", t2);
//
//  fs << "calibration_time" << buf;
//
//  if (!rvecs.empty() || !reprojErrs.empty())
//    fs << "nbFrames" << (int) std::max(rvecs.size(), reprojErrs.size());
//  fs << "image_width" << imageSize.width;
//  fs << "image_height" << imageSize.height;
//  fs << "board_width" << boardSize.width;
//  fs << "board_height" << boardSize.height;
//  fs << "square_size" << squareSize;
//
//  if (cvCalibFlags & CV_CALIB_FIX_ASPECT_RATIO)
//    fs << "aspectRatio" << aspectRatio;
//
//  if (cvCalibFlags != 0)
//  {
//    sprintf(buf, "flags: %s%s%s%s",
//            cvCalibFlags & CV_CALIB_USE_INTRINSIC_GUESS ? "+use_intrinsic_guess" : "",
//            cvCalibFlags & CV_CALIB_FIX_ASPECT_RATIO ? "+fix_aspectRatio" : "",
//            cvCalibFlags & CV_CALIB_FIX_PRINCIPAL_POINT ? "+fix_principal_point" : "",
//            cvCalibFlags & CV_CALIB_ZERO_TANGENT_DIST ? "+zero_tangent_dist" : "");
//    cvWriteComment(*fs, buf, 0);
//  }
//
//  fs << "flags" << cvCalibFlags;
//
//  fs << "camera_matrix" << cameraMatrix;
//  fs << "distortion_coefficients" << distCoeffs;
//
//  fs << "avg_reprojection_error" << totalAvgErr;
//  if (!reprojErrs.empty())
//    fs << "per_view_reprojection_errors" << cv::Mat(reprojErrs);
//
//  if (!rvecs.empty() && !tvecs.empty())
//  {
//    CV_Assert(rvecs[0].type() == tvecs[0].type());
//    cv::Mat bigmat((int) rvecs.size(), 6, rvecs[0].type());
//    for (int i = 0; i < (int) rvecs.size(); i++)
//    {
//      cv::Mat r = bigmat(cv::Range(i, i + 1), cv::Range(0, 3));
//      cv::Mat t = bigmat(cv::Range(i, i + 1), cv::Range(3, 6));
//
//      CV_Assert(rvecs[i].rows == 3 && rvecs[i].cols == 1);
//      CV_Assert(tvecs[i].rows == 3 && tvecs[i].cols == 1);
//      //*.t() is MatExpr (not Mat) so we can use assignment operator
//      r = rvecs[i].t();
//      t = tvecs[i].t();
//    }
//    cvWriteComment(*fs, "a set of 6-tuples (rotation vector + translation vector) for each view", 0);
//    fs << "extrinsic_parameters" << bigmat;
//  }
//
//  if (!imagePoints.empty())
//  {
//    cv::Mat imagePtMat((int) imagePoints.size(), (int) imagePoints[0].size(), CV_32FC2);
//    for (int i = 0; i < (int) imagePoints.size(); i++)
//    {
//      cv::Mat r = imagePtMat.row(i).reshape(2, imagePtMat.cols);
//      cv::Mat imgpti(imagePoints[i]);
//      imgpti.copyTo(r);
//    }
//    fs << "image_points" << imagePtMat;
//  }
//  const std::string txtfilename = filename.substr(0, filename.find_last_of(".")) + ".cal.txt";
//  saveCameraParamsToPlainTxt(txtfilename, imageSize, cameraMatrix, distCoeffs);
//}


///**
// * @brief This function computes cell indexes per image.
// *
// * @param[in] imagePoints
// * @param[out] cellIndexesPerImage
// * @param[in] imageSize
// * @param[in] calibGridSize
// */
//void precomputeCellIndexes(const std::vector<std::vector<cv::Point2f> >& imagePoints,
//                           std::vector<std::vector<std::size_t> >& cellIndexesPerImage,
//                           const cv::Size& imageSize, const std::size_t calibGridSize)
//{
//  float cellWidth = float(imageSize.width) / float(calibGridSize);
//  float cellHeight = float(imageSize.height) / float(calibGridSize);
//
//  for (const auto& pointbuf : imagePoints)
//  {
//    std::vector<std::size_t> imageCellIndexes;
//    // Points repartition in image
//    for (cv::Point2f point : pointbuf)
//    {
//      // Compute the index of the point
//      std::size_t cellPointX = std::floor(point.x / cellWidth);
//      std::size_t cellPointY = std::floor(point.y / cellHeight);
//      std::size_t cellIndex = cellPointY * calibGridSize + cellPointX;
//      imageCellIndexes.push_back(cellIndex);
//    }
//    cellIndexesPerImage.push_back(imageCellIndexes);
//  }
//}
//
///**
// * @brief This function counts the number of points in each cell of the grid.
// *
// * @param[in] imagesIndexes
// * @param[in] cellIndexesPerImage
// * @param[out] cellsWeight
// * @param[in] calibGridSize
// */
//void computeCellsWeight(const std::vector<std::size_t>& imagesIndexes,
//                        const std::vector<std::vector<std::size_t> >& cellIndexesPerImage,
//                        std::map<std::size_t, std::size_t>& cellsWeight, const std::size_t calibGridSize)
//{
//  //Init cell's weight to 0
//  for (std::size_t i = 0; i < calibGridSize * calibGridSize; ++i)
//    cellsWeight[i] = 0;
//
//  // Add weight into cells
//  for (std::size_t i = 0; i < imagesIndexes.size(); ++i)
//  {
//    std::vector<std::size_t> uniqueCellIndexes = cellIndexesPerImage[imagesIndexes[i]];
//    std::sort(uniqueCellIndexes.begin(), uniqueCellIndexes.end());
//    auto last = std::unique(uniqueCellIndexes.begin(), uniqueCellIndexes.end());
//    uniqueCellIndexes.erase(last, uniqueCellIndexes.end());
//
//    for (std::size_t cellIndex : uniqueCellIndexes)
//    {
//      ++cellsWeight[cellIndex];
//    }
//  }
//}
//
///**
// * @brief This function computes the score of each image.
// *
// * @param[out] imageScores
// * @param[in] remainingImagesIndexes
// * @param[in] cellIndexesPerImage
// * @param[in] cellsWeight
// */
//void computeImageScores(std::vector<std::pair<float, std::size_t> >& imageScores,
//                        const std::vector<std::size_t>& remainingImagesIndexes,
//                        const std::vector<std::vector<std::size_t> >& cellIndexesPerImage,
//                        const std::map<std::size_t, std::size_t>& cellsWeight)
//{
//  // Compute the score of each image
//  for (std::size_t i = 0; i < remainingImagesIndexes.size(); ++i)
//  {
//    const std::vector<std::size_t>& imageCellIndexes = cellIndexesPerImage[remainingImagesIndexes[i]];
//    float imageScore = 0;
//    for (std::size_t cellIndex : imageCellIndexes)
//    {
//      imageScore += cellsWeight.at(cellIndex);
//    }
//    // Normalize by the number of checker items.
//    // If the detector support occlusions of the checker the number of items may vary.
//    imageScore /= float(imageCellIndexes.size());
//    imageScores.emplace_back(imageScore, remainingImagesIndexes[i]);
//  }
//}

///**
// * @brief This function selects the best images based on repartition in images of the calibration landmarks.
// *
// * @param[in] imagePoints
// * @param[in] imageSize
// * @param[in,out] remainingImagesIndexes
// * @param[in] maxCalibFrames
// * @param[in] validFrames
// * @param[out] calibImageScore
// * @param[out] calibInputFrames
// * @param[out] calibImagePoints
// * @param[in] calibGridSize
// */
//void selectBestImages(const std::vector<std::vector<cv::Point2f> >& imagePoints,
//                      const cv::Size& imageSize,
//                      std::vector<std::size_t>& remainingImagesIndexes,
//                      const std::size_t& maxCalibFrames,
//                      const std::vector<std::size_t>& validFrames,
//                      std::vector<float>& calibImageScore,
//                      std::vector<std::size_t>& calibInputFrames,
//                      std::vector<std::vector<cv::Point2f> >& calibImagePoints,
//                      const std::size_t calibGridSize)
//{
//  std::vector<std::vector<std::size_t> > cellIndexesPerImage;
//
//  // Precompute cell indexes per image
//  precomputeCellIndexes(imagePoints, cellIndexesPerImage, imageSize, calibGridSize);
//
//  // Init with 0, 1, 2, ...
//  for (std::size_t i = 0; i < remainingImagesIndexes.size(); ++i)
//    remainingImagesIndexes[i] = i;
//
//  std::vector<std::size_t> bestImagesIndexes;
//  if (maxCalibFrames < validFrames.size())
//  {
//    while (bestImagesIndexes.size() < maxCalibFrames )
//    {
//      std::map<std::size_t, std::size_t> cellsWeight;
//      std::vector<std::pair<float, std::size_t> > imageScores;
//      // Count points in each cell of the grid
//      if (bestImagesIndexes.empty())
//        computeCellsWeight(remainingImagesIndexes, cellIndexesPerImage, cellsWeight, calibGridSize);
//      else
//        computeCellsWeight(bestImagesIndexes, cellIndexesPerImage, cellsWeight, calibGridSize);
//
//      computeImageScores(imageScores, remainingImagesIndexes, cellIndexesPerImage, cellsWeight);
//
//      // Find best score
//      std::size_t bestImageIndex = std::numeric_limits<std::size_t>::max();
//      float bestScore = std::numeric_limits<float>::max();
//      for (const auto& imageScore: imageScores)
//      {
//        if (imageScore.first < bestScore)
//        {
//          bestScore = imageScore.first;
//          bestImageIndex = imageScore.second;
//        }
//      }
//      auto eraseIt = std::find(remainingImagesIndexes.begin(), remainingImagesIndexes.end(), bestImageIndex);
//      assert(bestScore != std::numeric_limits<float>::max());
//      assert(eraseIt != remainingImagesIndexes.end());
//      remainingImagesIndexes.erase(eraseIt);
//      bestImagesIndexes.push_back(bestImageIndex);
//      calibImageScore.push_back(bestScore);
//    }
//  }
//  else
//  {
//    std::cout << "Info: Less valid frames (" << validFrames.size() << ") than specified maxCalibFrames (" << maxCalibFrames << ")." << std::endl;
//    bestImagesIndexes = validFrames;
//    std::map<std::size_t, std::size_t> cellsWeight;
//    computeCellsWeight(remainingImagesIndexes, cellIndexesPerImage, cellsWeight, calibGridSize);
//    std::vector<std::pair<float, std::size_t> > imageScores;
//    computeImageScores(imageScores, remainingImagesIndexes, cellIndexesPerImage, cellsWeight);
//    for(auto imgScore: imageScores)
//    {
//      calibImageScore.push_back(imgScore.first);
//    }
//  }
//
//  assert(bestImagesIndexes.size() == std::min(maxCalibFrames, validFrames.size()));
//
//  for(std::size_t i = 0; i < bestImagesIndexes.size(); ++i)
//  {
//    const std::size_t origI = bestImagesIndexes[i];
//    calibImagePoints.push_back(imagePoints[origI]);
//    calibInputFrames.push_back(validFrames[origI]);
//  }
//}

std::istream& operator>>(std::istream &in, openMVG::patternDetect::Pattern &pattern)
{
  std::string token;
  in >> token;
  boost::to_upper(token);

  if (token == "CHESSBOARD")
    pattern = openMVG::patternDetect::Pattern::CHESSBOARD;
  else if (token == "CIRCLES")
    pattern = openMVG::patternDetect::Pattern::CIRCLES_GRID;
  else if (token == "ASYMMETRIC_CIRCLES")
    pattern = openMVG::patternDetect::Pattern::ASYMMETRIC_CIRCLES_GRID;
#ifdef HAVE_CCTAG
    //  else if (token == "CCTAG")
    //    pattern = CCTAG_GRID;
#endif
  else
    throw boost::program_options::invalid_option_value(std::string("Invalid pattern: ") + token);
  return in;
}

int main(int argc, char** argv)
{
  // Command line arguments
  bfs::path inputPath;
  std::string outputFilename;
  std::string debugSelectedImgFolder;
  std::string debugRejectedImgFolder;
  std::vector<std::size_t> checkerboardSize;
  openMVG::patternDetect::Pattern pattern = openMVG::patternDetect::Pattern::CHESSBOARD;
  std::size_t maxNbFrames = 0;
  std::size_t maxCalibFrames = 100;
  std::size_t calibGridSize = 10;
  std::size_t nbDistortionCoef = 3;
  std::size_t minInputFrames = 10;
  double squareSize = 1.0;
  double maxTotalAvgErr = 0.1;

  std::clock_t startAlgo = std::clock();
  double durationAlgo;

  po::options_description desc("\n\nThis program is used to calibrate a camera from a dataset of images.\n");
  desc.add_options()
          ("help,h", "Produce help message.\n")
          ("input,i", po::value<bfs::path>(&inputPath)->required(),
           "Input images in one of the following form:\n"
           " - folder containing images\n"
           " - image sequence like /path/to/seq.@.jpg\n"
           " - video file\n")
          ("output,o", po::value<std::string>(&outputFilename)->required(),
           "Output filename for intrinsic [and extrinsic] parameters.\n")
//          ("pattern,p", po::value<openMVG::patternDetect::Pattern>(&pattern)->default_value(pattern),
//           "Type of pattern: 'chessboard', 'circles', 'asymmetric_circles'"
//#ifdef HAVE_CCTAG
//          " or 'cctag'"
//#endif
//          ".\n")
          ("size,s", po::value<std::vector < std::size_t >> (&checkerboardSize)->multitoken(),
           "Number of inner corners per one of board dimension like W H.\n")
          ("squareSize", po::value<double> (&squareSize)->default_value(squareSize),
           "Size of the grid's square cells (mm).\n")
          ("nbDistortionCoef,r", po::value<std::size_t>(&nbDistortionCoef)->default_value(nbDistortionCoef),
           "Number of distortion coefficient.\n")
          ("maxFrames", po::value<std::size_t>(&maxNbFrames)->default_value(maxNbFrames),
           "Maximal number of frames to extract from the video file.\n")
          ("maxCalibFrames", po::value<std::size_t>(&maxCalibFrames)->default_value(maxCalibFrames),
           "Maximal number of frames to use to calibrate from the selected frames.\n")
          ("calibGridSize", po::value<std::size_t>(&calibGridSize)->default_value(calibGridSize),
           "Define the number of cells per edge.\n")
          ("minInputFrames", po::value<std::size_t>(&minInputFrames)->default_value(minInputFrames),
           "Minimal number of frames to limit the refinement loop.\n")
          ("maxTotalAvgErr,e", po::value<double>(&maxTotalAvgErr)->default_value(maxTotalAvgErr),
           "Max Total Average Error.\n")
          ("debugRejectedImgFolder", po::value<std::string>(&debugRejectedImgFolder)->default_value(""),
           "Folder to export delete images during the refinement loop.\n")
          ("debugSelectedImgFolder,d", po::value<std::string>(&debugSelectedImgFolder)->default_value(""),
           "Folder to export debug images.\n")
          ;

  po::variables_map vm;
  int cvCalibFlags = 0;

  try
  {
    po::store(po::parse_command_line(argc, argv, desc), vm);

    if (vm.count("help") || (argc == 1))
    {
      std::cout << desc << std::endl;
      return EXIT_SUCCESS;
    }
    
    cvCalibFlags |= CV_CALIB_ZERO_TANGENT_DIST;
    if (nbDistortionCoef < 1 || nbDistortionCoef > 6)
      throw boost::program_options::invalid_option_value(std::string("Only supports 2 or 3 radial coefs: ") + std::to_string(nbDistortionCoef));
    const std::array<int, 6> fixDistortionCoefs = {CV_CALIB_FIX_K1, CV_CALIB_FIX_K2, CV_CALIB_FIX_K3, CV_CALIB_FIX_K4, CV_CALIB_FIX_K5, CV_CALIB_FIX_K6};
    for (int i = nbDistortionCoef; i < 6; ++i)
      cvCalibFlags |= fixDistortionCoefs[i];

    po::notify(vm);
  }
  catch (boost::program_options::required_option& e)
  {
    std::cerr << "ERROR: " << e.what() << std::endl;
    std::cout << "Usage:\n\n" << desc << std::endl;
    return EXIT_FAILURE;
  }
  catch (boost::program_options::error& e)
  {
    std::cerr << "ERROR: " << e.what() << std::endl;
    std::cout << "Usage:\n\n" << desc << std::endl;
    return EXIT_FAILURE;
  }

  if (checkerboardSize.size() != 2)
    throw std::logic_error("The size of the checkerboard is not defined");

  if (maxCalibFrames > maxNbFrames || minInputFrames > maxCalibFrames)
  {
    throw std::logic_error("Check the value for maxFrames, maxCalibFrames & minInputFrames. It must be decreasing.");
  }

  bool writeExtrinsics = false;
  bool writePoints = false;
  float aspectRatio = 1.f;
  cv::Mat cameraMatrix;
  cv::Mat distCoeffs;

  cv::Size boardSize(checkerboardSize[0], checkerboardSize[1]);
  cv::Size imageSize(0, 0);

  std::vector<std::vector<cv::Point2f> > imagePoints;

  std::clock_t start = std::clock();
  double duration;

  // create the feedProvider
  openMVG::dataio::FeedProvider feed(inputPath.string());
  if (!feed.isInit())
  {
    std::cerr << "ERROR while initializing the FeedProvider!" << std::endl;
    return EXIT_FAILURE;
  }
  openMVG::image::Image<unsigned char> imageGrey;
  openMVG::cameras::Pinhole_Intrinsic_Radial_K3 queryIntrinsics;
  bool hasIntrinsics = false;
  std::string currentImgName;
  std::size_t iInputFrame = 0;
  std::vector<std::size_t> validFrames;
  double step = 1.0;

  // Compute the discretization's step
  if (maxNbFrames)
  {
    if (feed.nbFrames() < maxNbFrames)
      step = feed.nbFrames() / (double) maxNbFrames;
  }

  while (feed.readImage(imageGrey, queryIntrinsics, currentImgName, hasIntrinsics) && iInputFrame < maxNbFrames)
  {
    std::size_t currentFrame = std::floor(iInputFrame * step);
    cv::Mat viewGray;
    cv::eigen2cv(imageGrey.GetMat(), viewGray);

    // Check image is correctly loaded
    if (viewGray.size() == cv::Size(0, 0))
    {
      throw std::runtime_error(std::string("Invalid image: ") + currentImgName);
    }
    // Check image size is always the same
    if (imageSize == cv::Size(0, 0))
    {
      // First image: initialize the image size.
      imageSize = viewGray.size();
    }
    // Check image resolutions are always the same
    else if (imageSize != viewGray.size())
    {
      throw std::runtime_error(std::string("You cannot mix multiple image resolutions during the camera calibration. See image file: ") + currentImgName);
    }

    std::vector<cv::Point2f> pointbuf;
    bool found;
    std::cout << "[" << currentFrame << "/" << maxNbFrames << "]" << std::endl;

    // Find the chosen pattern in images
    openMVG::patternDetect::findPattern(pattern, found, viewGray, boardSize, pointbuf);

    if (found)
    {
      validFrames.push_back(currentFrame);
      imagePoints.push_back(pointbuf);
    }

    ++iInputFrame;
    feed.goToFrame(std::floor(currentFrame));
  }

  duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
  std::cout << "find points duration: " << duration << std::endl;
  std::cout << "Grid detected in " << imagePoints.size() << " images on " << iInputFrame << " input images." << std::endl;

  if (imagePoints.empty())
    throw std::logic_error("No checkerboard detected.");

  std::vector<std::size_t> remainingImagesIndexes(validFrames.size());
  std::vector<float> calibImageScore;
  std::vector<std::size_t> calibInputFrames;
  std::vector<std::vector<cv::Point2f> > calibImagePoints;

  // Select best images based on repartition in images of the calibration landmarks
  openMVG::bestImages::selectBestImages(imagePoints, imageSize, remainingImagesIndexes, maxCalibFrames,
                   validFrames, calibImageScore, calibInputFrames, calibImagePoints, calibGridSize);

  std::vector<std::vector<cv::Point3f> > calibObjectPoints;
  start = std::clock();

  openMVG::patternDetect::computeObjectPoints(boardSize, pattern, squareSize, calibImagePoints, calibObjectPoints);

  double totalAvgErr = 0;
  std::vector<cv::Mat> rvecs;
  std::vector<cv::Mat> tvecs;
  std::vector<float> reprojErrs;
  std::vector<std::size_t> rejectInputFrames;
  
  // Refinement loop of the calibration
  openMVG::calibration::calibrationIterativeOptimization(calibImagePoints, calibObjectPoints, imageSize, aspectRatio,
                        cvCalibFlags, cameraMatrix, distCoeffs, rvecs, tvecs, reprojErrs,
                        totalAvgErr, maxTotalAvgErr, minInputFrames, calibInputFrames,
                        calibImageScore, rejectInputFrames);

  duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
  std::cout << "Calibration duration: " << duration << std::endl;

  openMVG::exportData::saveCameraParams(outputFilename, imageSize,
                   boardSize, squareSize, aspectRatio,
                   cvCalibFlags, cameraMatrix, distCoeffs,
                   writeExtrinsics ? rvecs : std::vector<cv::Mat>(),
                   writeExtrinsics ? tvecs : std::vector<cv::Mat>(),
                   writeExtrinsics ? reprojErrs : std::vector<float>(),
                   writePoints ? calibImagePoints : std::vector<std::vector<cv::Point2f> >(),
                   totalAvgErr);

  openMVG::exportData::exportDebug(debugSelectedImgFolder, debugRejectedImgFolder,
              feed, calibInputFrames, rejectInputFrames, remainingImagesIndexes,
              cameraMatrix, distCoeffs, imageSize);

  durationAlgo = (std::clock() - startAlgo) / (double) CLOCKS_PER_SEC;
  std::cout << "Total duration: " << durationAlgo << std::endl;

  return 0;
}