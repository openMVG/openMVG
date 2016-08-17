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

int main(int argc, char** argv)
{
  // Command line arguments
  bfs::path inputPath;
  std::string outputFilename;
  std::string debugSelectedImgFolder;
  std::string debugRejectedImgFolder;
  std::vector<std::size_t> checkerboardSize;
  openMVG::calibration::Pattern pattern = openMVG::calibration::Pattern::CHESSBOARD;
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
          ("pattern,p", po::value<openMVG::calibration::Pattern>(&pattern)->default_value(pattern),
           "Type of pattern: 'chessboard', 'circles', 'asymmetric_circles'"
            #ifdef HAVE_CCTAG
//                      " or 'cctag'"
            #endif
          ".\n")
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
    std::cout << "[" << currentFrame << "/" << maxNbFrames << "]" << std::endl;

    // Find the chosen pattern in images
    const bool found = openMVG::calibration::findPattern(pattern, viewGray, boardSize, pointbuf);

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
  openMVG::calibration::selectBestImages(imagePoints, imageSize, remainingImagesIndexes, maxCalibFrames,
                                         validFrames, calibImageScore, calibInputFrames, calibImagePoints, calibGridSize);

  std::vector<std::vector<cv::Point3f> > calibObjectPoints;
  start = std::clock();

  openMVG::calibration::computeObjectPoints(boardSize, pattern, squareSize, calibImagePoints, calibObjectPoints);

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

  openMVG::calibration::saveCameraParams(outputFilename, imageSize,
                                         boardSize, squareSize, aspectRatio,
                                         cvCalibFlags, cameraMatrix, distCoeffs,
                                         writeExtrinsics ? rvecs : std::vector<cv::Mat>(),
                                         writeExtrinsics ? tvecs : std::vector<cv::Mat>(),
                                         writeExtrinsics ? reprojErrs : std::vector<float>(),
                                         writePoints ? calibImagePoints : std::vector<std::vector<cv::Point2f> >(),
                                         totalAvgErr);

  openMVG::calibration::exportDebug(debugSelectedImgFolder, debugRejectedImgFolder,
                                    feed, calibInputFrames, rejectInputFrames, remainingImagesIndexes,
                                    cameraMatrix, distCoeffs, imageSize);

  durationAlgo = (std::clock() - startAlgo) / (double) CLOCKS_PER_SEC;
  std::cout << "Total duration: " << durationAlgo << std::endl;

  return 0;
}
