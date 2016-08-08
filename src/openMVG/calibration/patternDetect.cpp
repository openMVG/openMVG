#include "patternDetect.hpp"

#include <boost/lexical_cast.hpp>
#include <ctime>
#include <cctype>
#include <stdexcept>
#include <iostream>

namespace openMVG {
namespace patternDetect {

int findPattern(const Pattern& pattern, bool& found, const cv::Mat& viewGray, const cv::Size& boardSize, std::vector<cv::Point2f>& pointbuf)
{
  std::clock_t startCh;
  double durationCh;

  switch (pattern)
  {
  case CHESSBOARD:
    startCh = std::clock();

    found = cv::findChessboardCorners(viewGray, boardSize, pointbuf,
                                      CV_CALIB_CB_ADAPTIVE_THRESH | CV_CALIB_CB_FAST_CHECK | CV_CALIB_CB_NORMALIZE_IMAGE);
    durationCh = (std::clock() - startCh) / (double) CLOCKS_PER_SEC;
    std::cout << "Find chessboard corners' duration: " << durationCh << std::endl;
    startCh = std::clock();

    // improve the found corners' coordinate accuracy
    if (found)
      cv::cornerSubPix(viewGray, pointbuf, cv::Size(11, 11), cv::Size(-1, -1),
                       cv::TermCriteria(CV_TERMCRIT_EPS + CV_TERMCRIT_ITER, 30, 0.1));

    durationCh = (std::clock() - startCh) / (double) CLOCKS_PER_SEC;
    std::cout << "Refine chessboard corners' duration: " << durationCh << std::endl;
    break;

  case CIRCLES_GRID:
    startCh = std::clock();

    found = cv::findCirclesGrid(viewGray, boardSize, pointbuf);

    durationCh = (std::clock() - startCh) / (double) CLOCKS_PER_SEC;
    std::cout << "Find circles grid duration: " << durationCh << std::endl;
    break;

  case ASYMMETRIC_CIRCLES_GRID:
    startCh = std::clock();

    found = cv::findCirclesGrid(viewGray, boardSize, pointbuf, cv::CALIB_CB_ASYMMETRIC_GRID);

    durationCh = (std::clock() - startCh) / (double) CLOCKS_PER_SEC;
    std::cout << "Find asymmetric circles grid duration: " << durationCh << std::endl;
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
  return EXIT_SUCCESS;
}

void calcChessboardCorners(cv::Size boardSize, const float& squareSize,
                                  std::vector<cv::Point3f>& corners, Pattern pattern = Pattern::CHESSBOARD)
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

void computeObjectPoints(cv::Size boardSize,
                                Pattern pattern,
                                const float& squareSize,
                                const std::vector<std::vector<cv::Point2f> >& imagePoints,
                                std::vector<std::vector<cv::Point3f> >& objectPoints)
{
  std::vector<cv::Point3f> templateObjectPoints;

  // Generate the object points coordinates
  calcChessboardCorners(boardSize, squareSize, templateObjectPoints, pattern);

  // Assign the template to all items
  objectPoints.resize(imagePoints.size(), templateObjectPoints);
}

}//namespace patternDetect
}//namespace openMVG



