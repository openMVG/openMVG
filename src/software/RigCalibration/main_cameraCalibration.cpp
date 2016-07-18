#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/algorithm/string/case_conv.hpp>

#include <openMVG/dataio/FeedProvider.hpp>
#include <openMVG/cameras/Camera_undistort_image.hpp>
#include <openMVG/cameras/Camera_Pinhole_Radial.hpp>
#include <openMVG/image/image_io.hpp>

#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/core/eigen.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/calib3d/calib3d.hpp>

#include <stdio.h>
#include <time.h>
#include <ctime>
#include <cstdio>
#include <string>
#include <cctype>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <exception>

namespace bfs = boost::filesystem;
namespace po = boost::program_options;

enum Pattern
{
  CHESSBOARD,
  CIRCLES_GRID,
  ASYMMETRIC_CIRCLES_GRID
#ifdef HAVE_CCTAG
  , CCTAG_GRID
#endif
};

void exportDebug(openMVG::dataio::FeedProvider& feed, const std::string& debugFolder, 
                const std::vector<int>& validFrames, const cv::Mat& cameraMatrix, 
                const cv::Mat& distCoeffs, const cv::Size& imageSize)
{
  std::vector<int> export_params;
  openMVG::image::Image<unsigned char> inputImage;
  openMVG::image::Image<unsigned char> outputImage;
  std::string currentImgName;
  openMVG::cameras::Pinhole_Intrinsic_Radial_K3 queryIntrinsics;
  bool hasIntrinsics;
  cv::Mat rview;

  export_params.push_back(CV_IMWRITE_JPEG_QUALITY);
  export_params.push_back(100);

  std::cout << "Exporting undistorted images ..." << std::endl;
  for(int currentFrame : validFrames)
  {
    feed.goToFrame(currentFrame);
    feed.readImage(inputImage, queryIntrinsics, currentImgName, hasIntrinsics);

    // drawChessboardCorners(view, boardSize, cv::Mat(pointbuf), found);

    openMVG::cameras::Pinhole_Intrinsic_Radial_K3 camera(
            imageSize.width, imageSize.height,
            cameraMatrix.at<double>(0,0), cameraMatrix.at<double>(0,2), cameraMatrix.at<double>(1,2),
            distCoeffs.at<double>(0), distCoeffs.at<double>(1), distCoeffs.at<double>(2));

    openMVG::cameras::UndistortImage(inputImage, &camera, outputImage, openMVG::image::BLACK);
    const bfs::path imagePath = bfs::path(debugFolder) / (std::to_string(currentFrame) + "_undistort.png");
    const bool exportStatus = openMVG::image::WriteImage(imagePath.c_str(), outputImage);
    if(!exportStatus)
    {
      std::cerr << "Failed to export: " << imagePath << std::endl;
    }
  }
  std::cout << "... finished" << std::endl;
}

static double computeReprojectionErrors(
                                        const std::vector<std::vector<cv::Point3f> >& objectPoints,
                                        const std::vector<std::vector<cv::Point2f> >& imagePoints,
                                        const std::vector<cv::Mat>& rvecs, const std::vector<cv::Mat>& tvecs,
                                        const cv::Mat& cameraMatrix, const cv::Mat& distCoeffs,
                                        std::vector<float>& perViewErrors)
{
  std::vector<cv::Point2f> imagePoints2;
  int i, totalPoints = 0;
  double err, totalErr = 0;
  perViewErrors.resize(objectPoints.size());

  for (i = 0; i < (int) objectPoints.size(); i++)
  {
    cv::projectPoints(cv::Mat(objectPoints[i]), rvecs[i], tvecs[i],
                  cameraMatrix, distCoeffs, imagePoints2);
    err = cv::norm(cv::Mat(imagePoints[i]), cv::Mat(imagePoints2), CV_L2);
    int n = (int) objectPoints[i].size();
    perViewErrors[i] = (float) std::sqrt(err * err / n);
    totalErr += err*err;
    totalPoints += n;
  }

  return std::sqrt(totalErr / totalPoints);
}

static void calcChessboardCorners(cv::Size boardSize, float squareSize, 
                                  std::vector<cv::Point3f>& corners, Pattern pattern = CHESSBOARD)
{
  corners.resize(0);

  switch (pattern)
  {
  case CHESSBOARD:
  case CIRCLES_GRID:
    for (int i = 0; i < boardSize.height; i++)
      for (int j = 0; j < boardSize.width; j++)
        corners.push_back(cv::Point3f(float(j * squareSize),
                                  float(i * squareSize), 0));
    break;

  case ASYMMETRIC_CIRCLES_GRID:
    for (int i = 0; i < boardSize.height; i++)
      for (int j = 0; j < boardSize.width; j++)
        corners.push_back(cv::Point3f(float((2 * j + i % 2) * squareSize),
                                  float(i * squareSize), 0));
    break;

#ifdef HAVE_CCTAG
  case CCTAG_GRID:
    throw std::invalid_argument("CCTag not implemented.");
    break;
#endif

  default:
    throw std::invalid_argument("Unknown pattern type.");
  }
}

static bool runCalibration(std::vector<std::vector<cv::Point2f> > imagePoints,
                           cv::Size imageSize, cv::Size boardSize, Pattern pattern,
                           float squareSize, float aspectRatio,
                           int flags, cv::Mat& cameraMatrix, cv::Mat& distCoeffs,
                           std::vector<cv::Mat>& rvecs, std::vector<cv::Mat>& tvecs,
                           std::vector<float>& reprojErrs,
                           double& totalAvgErr)
{
  cameraMatrix = cv::Mat::eye(3, 3, CV_64F);
  if (flags & CV_CALIB_FIX_ASPECT_RATIO)
    cameraMatrix.at<double>(0, 0) = aspectRatio;

  distCoeffs = cv::Mat::zeros(8, 1, CV_64F);

  std::vector<std::vector<cv::Point3f> > objectPoints(1);
  
  std::clock_t startrC = std::clock();
  double durationrC;
  
  calcChessboardCorners(boardSize, squareSize, objectPoints[0], pattern);
  
  durationrC = ( std::clock() - startrC ) / (double) CLOCKS_PER_SEC;
  std::cout <<"  calcChessboardCorners duration: "<< durationrC << std::endl;
  
  objectPoints.resize(imagePoints.size(), objectPoints[0]);
 
  startrC = std::clock();
  
  double rms = cv::calibrateCamera(objectPoints, imagePoints, imageSize, cameraMatrix,
                               distCoeffs, rvecs, tvecs, flags | CV_CALIB_FIX_K4 | CV_CALIB_FIX_K5 | CV_CALIB_FIX_K6);
  durationrC = ( std::clock() - startrC ) / (double) CLOCKS_PER_SEC;
  std::cout <<"  calibrateCamera duration: "<< durationrC << std::endl;
  
  printf("RMS error reported by calibrateCamera: %g\n", rms);
  bool ok = cv::checkRange(cameraMatrix) && cv::checkRange(distCoeffs);

  startrC = std::clock();
  
  totalAvgErr = computeReprojectionErrors(objectPoints, imagePoints,
                                          rvecs, tvecs, cameraMatrix, distCoeffs, reprojErrs);
  durationrC = ( std::clock() - startrC ) / (double) CLOCKS_PER_SEC;
  std::cout <<"  computeReprojectionErrors duration: "<< durationrC << std::endl;

  return ok;
}

static void saveCameraParamsToPlainTxt(const std::string &filename,
                                       const cv::Size &imageSize,
                                       const cv::Mat& cameraMatrix,
                                       const cv::Mat& distCoeffs)
{
  std::ofstream fs(filename, std::ios::out);
  if (!fs.is_open())
  {
    std::cerr << "Unable to create the calibration file " << filename << std::endl;
    throw std::invalid_argument("Unable to create the calibration file " + filename);
  }

  // the structure of the file is
  // int #image width
  // int #image height
  // double #focal
  // double #ppx principal point x-coord
  // double #ppy principal point y-coord
  // double #k0
  // double #k1
  // double #k2
  fs << imageSize.width << std::endl;
  fs << imageSize.height << std::endl;
  if (cameraMatrix.type() == cv::DataType<double>::type)
  {
    fs << (cameraMatrix.at<double>(0, 0) + cameraMatrix.at<double>(1, 1)) / 2 << std::endl;
    fs << cameraMatrix.at<double>(0, 2) << std::endl;
    fs << cameraMatrix.at<double>(1, 2) << std::endl;
  }
  else
  {
    fs << (cameraMatrix.at<float>(0, 0) + cameraMatrix.at<float>(1, 1)) / 2 << std::endl;
    fs << cameraMatrix.at<float>(0, 2) << std::endl;
    fs << cameraMatrix.at<float>(1, 2) << std::endl;
  }
  if (distCoeffs.type() == cv::DataType<double>::type)
  {
    fs << distCoeffs.at<double>(0) << std::endl;
    fs << distCoeffs.at<double>(1) << std::endl;
    fs << distCoeffs.at<double>(2) << std::endl;
  }
  else
  {
    fs << distCoeffs.at<float>(0) << std::endl;
    fs << distCoeffs.at<float>(1) << std::endl;
    fs << distCoeffs.at<float>(2) << std::endl;
  }
  fs.close();
}

static void saveCameraParams(const std::string& filename,
                             cv::Size imageSize, cv::Size boardSize,
                             float squareSize, float aspectRatio, int flags,
                             const cv::Mat& cameraMatrix, const cv::Mat& distCoeffs,
                             const std::vector<cv::Mat>& rvecs, const std::vector<cv::Mat>& tvecs,
                             const std::vector<float>& reprojErrs,
                             const std::vector<std::vector<cv::Point2f> >& imagePoints,
                             double totalAvgErr)
{
  cv::FileStorage fs(filename, cv::FileStorage::WRITE);

  time_t tt;
  time(&tt);
  struct tm *t2 = localtime(&tt);
  char buf[1024];
  strftime(buf, sizeof (buf) - 1, "%c", t2);

  fs << "calibration_time" << buf;

  if (!rvecs.empty() || !reprojErrs.empty())
    fs << "nbFrames" << (int) std::max(rvecs.size(), reprojErrs.size());
  fs << "image_width" << imageSize.width;
  fs << "image_height" << imageSize.height;
  fs << "board_width" << boardSize.width;
  fs << "board_height" << boardSize.height;
  fs << "square_size" << squareSize;

  if (flags & CV_CALIB_FIX_ASPECT_RATIO)
    fs << "aspectRatio" << aspectRatio;

  if (flags != 0)
  {
    sprintf(buf, "flags: %s%s%s%s",
            flags & CV_CALIB_USE_INTRINSIC_GUESS ? "+use_intrinsic_guess" : "",
            flags & CV_CALIB_FIX_ASPECT_RATIO ? "+fix_aspectRatio" : "",
            flags & CV_CALIB_FIX_PRINCIPAL_POINT ? "+fix_principal_point" : "",
            flags & CV_CALIB_ZERO_TANGENT_DIST ? "+zero_tangent_dist" : "");
    cvWriteComment(*fs, buf, 0);
  }

  fs << "flags" << flags;

  fs << "camera_matrix" << cameraMatrix;
  fs << "distortion_coefficients" << distCoeffs;

  fs << "avg_reprojection_error" << totalAvgErr;
  if (!reprojErrs.empty())
    fs << "per_view_reprojection_errors" << cv::Mat(reprojErrs);

  if (!rvecs.empty() && !tvecs.empty())
  {
    CV_Assert(rvecs[0].type() == tvecs[0].type());
    cv::Mat bigmat((int) rvecs.size(), 6, rvecs[0].type());
    for (int i = 0; i < (int) rvecs.size(); i++)
    {
      cv::Mat r = bigmat(cv::Range(i, i + 1), cv::Range(0, 3));
      cv::Mat t = bigmat(cv::Range(i, i + 1), cv::Range(3, 6));

      CV_Assert(rvecs[i].rows == 3 && rvecs[i].cols == 1);
      CV_Assert(tvecs[i].rows == 3 && tvecs[i].cols == 1);
      //*.t() is MatExpr (not Mat) so we can use assignment operator
      r = rvecs[i].t();
      t = tvecs[i].t();
    }
    cvWriteComment(*fs, "a set of 6-tuples (rotation vector + translation vector) for each view", 0);
    fs << "extrinsic_parameters" << bigmat;
  }

  if (!imagePoints.empty())
  {
    cv::Mat imagePtMat((int) imagePoints.size(), (int) imagePoints[0].size(), CV_32FC2);
    for (int i = 0; i < (int) imagePoints.size(); i++)
    {
      cv::Mat r = imagePtMat.row(i).reshape(2, imagePtMat.cols);
      cv::Mat imgpti(imagePoints[i]);
      imgpti.copyTo(r);
    }
    fs << "image_points" << imagePtMat;
  }
  const std::string txtfilename = filename.substr(0, filename.find_last_of(".")) + ".cal.txt";
  saveCameraParamsToPlainTxt(txtfilename, imageSize, cameraMatrix, distCoeffs);
}

std::istream& operator>> (std::istream &in, Pattern &pattern)
{
  std::string token;
  in >> token;
  boost::to_upper(token);

  if (token == "CHESSBOARD")
    pattern = CHESSBOARD;
  else if (token == "CIRCLES")
    pattern = CIRCLES_GRID;
  else if (token == "ASYMMETRIC_CIRCLES")
    pattern = ASYMMETRIC_CIRCLES_GRID;
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
  bfs::path inputPath;
  std::string outputFilename;
  std::string debugFolder;
  std::vector<std::size_t> checkerboardSize;
  Pattern pattern;
  unsigned int maxNbFrames = 0;
  unsigned int nbRadialCoef = 3;
  bool writeExtrinsics = false;
  bool writePoints = false;
  int flags = 0;
  float squareSize = 1.f;
  float aspectRatio = 1.f;
  cv::Mat cameraMatrix;
  cv::Mat distCoeffs;
  
  std::clock_t startAlgo = std::clock();
  double durationAlgo;
  
  po::options_description desc("This program is used to calibrate a camera from a dataset of images.\n");
  desc.add_options()
          ("help,h", "Produce help message.\n")
          ("input,i", po::value<bfs::path>(&inputPath)->required(), 
                      "Input images in one of the following form:\n"
                      " - folder containing images\n"
                      " - image sequence like /path/to/seq.@.jpg\n"
                      " - video file\n")
          ("output,o", po::value<std::string>(&outputFilename)->required(), 
                      "Output filename for intrinsic [and extrinsic] parameters.\n")
          ("pattern,p", po::value<Pattern>(&pattern)->default_value(CHESSBOARD),
                      "Type of pattern: 'chessboard', 'circles', 'asymmetric_circles'"
                      #ifdef HAVE_CCTAG
                        " or 'cctag'"
                      #endif
                      ".\n")
          ("size,s", po::value<std::vector<std::size_t>>(&checkerboardSize)->multitoken(),
                      "Number of inner corners per one of board dimension like W H.\n")
          ("maxFrames,m", po::value<unsigned int>(&maxNbFrames)->default_value(0),
                      "Number of maximal frames to use for calibration from a video file.\n")
          ("nRadialCoef,r", po::value<unsigned int>(&nbRadialCoef)->default_value(3), 
                      "Number of radial distortion coefficient.\n")
          ("debugFolder,d", po::value<std::string>(&debugFolder)->default_value(""),
                      "Folder to export debug images.\n")
  ;
  
  po::variables_map vm;
  
  try
  {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    
    if (vm.count("help") || (argc == 1))
    {
      std::cout << desc << std::endl;
      return EXIT_SUCCESS;
    }
    
    flags |= CV_CALIB_ZERO_TANGENT_DIST;
    if(nbRadialCoef < 1 || nbRadialCoef > 6)
      throw boost::program_options::invalid_option_value(std::string("Only supports 2 or 3 radial coefs: ") + std::to_string(nbRadialCoef));
    const std::array<int, 6> fixRadialCoefs = {CV_CALIB_FIX_K1, CV_CALIB_FIX_K2, CV_CALIB_FIX_K3, CV_CALIB_FIX_K4, CV_CALIB_FIX_K5, CV_CALIB_FIX_K6};
    for(int i = nbRadialCoef; i < 6; ++i)
      flags |= fixRadialCoefs[i];
    
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

  if(checkerboardSize.size() != 2)
    throw std::logic_error("The size of the checkerboard is not defined;");;
  cv::Size boardSize(checkerboardSize[0], checkerboardSize[1]);
  cv::Size imageSize(0, 0);

  std::vector<std::vector<cv::Point2f> > imagePoints;

  std::clock_t start = std::clock();
  double duration;

  // create the feedProvider
  openMVG::dataio::FeedProvider feed(inputPath.string());
  if(!feed.isInit())
  {
    std::cerr << "ERROR while initializing the FeedProvider!" << std::endl;
    return EXIT_FAILURE;
  }
  openMVG::image::Image<unsigned char> imageGrey;
  openMVG::cameras::Pinhole_Intrinsic_Radial_K3 queryIntrinsics;
  bool hasIntrinsics = false;
  std::string currentImgName;
  int step = 1;
  if(maxNbFrames)
    step = std::floor(feed.nbFrames() / (double)maxNbFrames);
  unsigned int iInputFrame = 0;
  std::vector<int> validFrames;

  while(feed.readImage(imageGrey, queryIntrinsics, currentImgName, hasIntrinsics))
  {
    cv::Mat viewGray;
    cv::eigen2cv(imageGrey.GetMat(), viewGray);

    // Check image is correctly loaded
    if (viewGray.size() == cv::Size(0,0))
    {
      throw std::runtime_error(std::string("Invalid image: ") + currentImgName);
    }
    // Check image size is always the same
    if (imageSize == cv::Size(0,0))
    {
      // First image: initialize the image size.
      imageSize = viewGray.size();
    }
    else if (imageSize != viewGray.size())
    {
      throw std::runtime_error(std::string("You cannot mix multiple image resolutions during the camera calibration. See image file: ") + currentImgName);
    }

    std::vector<cv::Point2f> pointbuf;

    std::clock_t startCh;
    double durationCh;
    bool found;
    switch (pattern)
    {
      case CHESSBOARD:
        startCh = std::clock();
        
        found = cv::findChessboardCorners(viewGray, boardSize, pointbuf,
                                      CV_CALIB_CB_ADAPTIVE_THRESH | CV_CALIB_CB_FAST_CHECK | CV_CALIB_CB_NORMALIZE_IMAGE);
        durationCh = ( std::clock() - startCh ) / (double) CLOCKS_PER_SEC;
        std::cout<< "find chessboard corners' duration: "<< durationCh << std::endl;
        startCh = std::clock();
        
        // improve the found corners' coordinate accuracy
        if (found)
          cv::cornerSubPix(viewGray, pointbuf, cv::Size(11, 11), cv::Size(-1, -1),
                       cv::TermCriteria(CV_TERMCRIT_EPS + CV_TERMCRIT_ITER, 30, 0.1));
        
        durationCh = ( std::clock() - startCh ) / (double) CLOCKS_PER_SEC;
        std::cout << "refine chessboard corners' duration: "<< durationCh << std::endl;
        break;
        
      case CIRCLES_GRID:
        startCh = std::clock();
        
        found = cv::findCirclesGrid(viewGray, boardSize, pointbuf);
        
        durationCh = ( std::clock() - startCh ) / (double) CLOCKS_PER_SEC;
        std::cout << "find circles grid duration: "<< durationCh << std::endl;
        break;
        
      case ASYMMETRIC_CIRCLES_GRID:
        startCh = std::clock();
        
        found = cv::findCirclesGrid(viewGray, boardSize, pointbuf, cv::CALIB_CB_ASYMMETRIC_GRID);
        
        durationCh = ( std::clock() - startCh ) / (double) CLOCKS_PER_SEC;
        std::cout << "find asymmetric circles grid duration: "<< durationCh << std::endl;
        break;

#ifdef HAVE_CCTAG
      case CCTAG_GRID:
        throw std::invalid_argument("CCTag calibration not implemented.");
        break;
#endif

      default:
        std::cerr << "Unknown pattern type" << std::endl;
        return EXIT_FAILURE;
    }

    if (found)
    {
      validFrames.push_back(iInputFrame * step);
      imagePoints.push_back(pointbuf);
    }
    
    ++iInputFrame;
    feed.goToFrame(iInputFrame * step);
  }
  
  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  std::cout << "find points duration: " << duration << std::endl;
  std::cout << "Grid detected in " << imagePoints.size() << " images on " << iInputFrame << " input images." << std::endl;

  if (imagePoints.empty())
    throw std::logic_error("No checkerboard detected.");

  std::vector<cv::Mat> rvecs;
  std::vector<cv::Mat> tvecs;
  std::vector<float> reprojErrs;
  double totalAvgErr = 0;

  start = std::clock();
  
  // TODO: Refinement Loop 
  bool calibSucceeded = runCalibration(imagePoints, imageSize, boardSize, pattern, squareSize,
                                       aspectRatio, flags, cameraMatrix, distCoeffs,
                                       rvecs, tvecs, reprojErrs, totalAvgErr);

  std::cout << " avg reprojection error = " << totalAvgErr << std::endl;
  std::cout << (calibSucceeded ? "Calibration succeeded" : "Calibration failed") << std::endl;
  
  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  std::cout << "Calibration duration: "<< duration << std::endl;

  if (!calibSucceeded)
    return -1;
  
  start = std::clock();
  
  saveCameraParams(outputFilename, imageSize,
                   boardSize, squareSize, aspectRatio,
                   flags, cameraMatrix, distCoeffs,
                   writeExtrinsics ? rvecs : std::vector<cv::Mat>(),
                   writeExtrinsics ? tvecs : std::vector<cv::Mat>(),
                   writeExtrinsics ? reprojErrs : std::vector<float>(),
                   writePoints ? imagePoints : std::vector<std::vector<cv::Point2f> >(),
                   totalAvgErr);
  
  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  std::cout << "saveCameraParams duration: "<< duration << std::endl;
  
  start = std::clock();
  if (!debugFolder.empty())
  {
    exportDebug(feed, debugFolder, validFrames, 
                cameraMatrix, distCoeffs, imageSize);
  }
  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  std::cout << "exportDebug duration: "<< duration << std::endl;
  
  durationAlgo = ( std::clock() - startAlgo ) / (double) CLOCKS_PER_SEC;
  std::cout << "total duration: "<< durationAlgo << std::endl;
  
  return 0;
}