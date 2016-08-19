#include "patternDetect.hpp"

#include <openMVG/system/timer.hpp>

#include <boost/program_options.hpp>
#include <boost/algorithm/string/case_conv.hpp>

#include <boost/lexical_cast.hpp>

#include <string>
#include <ctime>
#include <cctype>
#include <stdexcept>
#include <iostream>

namespace openMVG{
namespace calibration{

std::istream& operator>>(std::istream &stream, Pattern &pattern)
{
  std::string token;
  stream >> token;
  boost::to_upper(token);

  if (token == "CHESSBOARD")
    pattern = openMVG::calibration::Pattern::CHESSBOARD;
  else if (token == "CIRCLES")
    pattern = openMVG::calibration::Pattern::CIRCLES_GRID;
  else if (token == "ASYMMETRIC_CIRCLES")
    pattern = openMVG::calibration::Pattern::ASYMMETRIC_CIRCLES_GRID;
#ifdef HAVE_CCTAG
  else if (token == "CCTAG")
    pattern = CCTAG_GRID;
#endif
  else
    throw boost::program_options::invalid_option_value(std::string("Invalid pattern: ") + token);
  return stream;
}

std::ostream& operator<<(std::ostream &stream, const Pattern pattern)
{
  switch(pattern)
  {
    case openMVG::calibration::Pattern::CHESSBOARD:
      stream << "CHESSBOARD";
      break;
    case openMVG::calibration::Pattern::CIRCLES_GRID:
      stream << "CIRCLES_GRID";
      break;
    case openMVG::calibration::Pattern::ASYMMETRIC_CIRCLES_GRID:
      stream << "ASYMMETRIC_CIRCLES_GRID";
      break;
#ifdef HAVE_CCTAG
    case openMVG::calibration::Pattern::CCTAG_GRID:
      stream << "CCTAG_GRID";
      break;
#endif
  }
  return stream;
}

bool findPattern(const Pattern& pattern, const cv::Mat& viewGray, const cv::Size& boardSize, std::vector<cv::Point2f>& pointbuf)
{
  bool found = false;
  system::Timer durationCh;

  switch (pattern)
  {
  case CHESSBOARD:

    found = cv::findChessboardCorners(viewGray, boardSize, pointbuf,
                                      CV_CALIB_CB_ADAPTIVE_THRESH | CV_CALIB_CB_FAST_CHECK | CV_CALIB_CB_NORMALIZE_IMAGE);
    std::cout << "Find chessboard corners' duration: " << durationCh.elapsedMs() << "ms" << std::endl;

    // improve the found corners' coordinate accuracy
    if (found)
    {
      durationCh.reset();
      cv::cornerSubPix(viewGray, pointbuf, cv::Size(11, 11), cv::Size(-1, -1),
                       cv::TermCriteria(CV_TERMCRIT_EPS + CV_TERMCRIT_ITER, 30, 0.1));

      std::cout << "Refine chessboard corners' duration: " << durationCh.elapsedMs() << "ms" << std::endl;
    }
    break;

  case CIRCLES_GRID:

    found = cv::findCirclesGrid(viewGray, boardSize, pointbuf);

    std::cout << "Find circles grid duration: " << durationCh.elapsedMs() << "ms" << std::endl;
    break;

  case ASYMMETRIC_CIRCLES_GRID:

    found = cv::findCirclesGrid(viewGray, boardSize, pointbuf, cv::CALIB_CB_ASYMMETRIC_GRID);

    std::cout << "Find asymmetric circles grid duration: " << durationCh.elapsedMs() << "ms" << std::endl;
    break;

#ifdef HAVE_CCTAG
  case CCTAG_GRID:
    throw std::invalid_argument("CCTag calibration not implemented.");
    break;
#endif

  default:
    throw std::logic_error("LensCalibration: Unknown pattern type.");
  }
  return found;
}

void calcChessboardCorners(const cv::Size& boardSize, const float& squareSize,
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

void computeObjectPoints(const cv::Size& boardSize,
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

}//namespace calibration
}//namespace openMVG



