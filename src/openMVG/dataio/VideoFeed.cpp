/* 
 * File:   VideoFeed.cpp
 * Author: sgaspari
 * 
 * Created on September 28, 2015, 10:35 AM
 */

#if HAVE_OPENCV

#include "VideoFeed.hpp"

#include <openMVG/image/image_converter.hpp>

#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/core/eigen.hpp>
#include <opencv2/highgui.hpp>

#include <iostream>
#include <exception>

namespace openMVG{
namespace dataio{

class VideoFeed::FeederImpl
{
public:
  FeederImpl() : _isInit(false) { }
  
  FeederImpl(const std::string &videoPath, const std::string &calibPath);
  
  bool isInit() const {return _isInit;}
  
  bool readImage(image::Image<unsigned char> &imageGray,
                     cameras::Pinhole_Intrinsic_Radial_K3 &camIntrinsics,
                     std::string &mediaPath,
                     bool &hasIntrinsics);
  
  bool goToFrame(const unsigned int frame);
  
  bool goToNextFrame();
  
  std::size_t nbFrames() const;
  
private:
  bool _isInit;
  bool _withIntrinsics;
  std::string _videoPath;
  cv::VideoCapture _videoCapture;
  cameras::Pinhole_Intrinsic_Radial_K3 _camIntrinsics;
};


VideoFeed::FeederImpl::FeederImpl(const std::string &videoPath, const std::string &calibPath)
: _isInit(false), _withIntrinsics(false), _videoPath(videoPath)
{
    // load the video
  _videoCapture.open(videoPath);
  if (!_videoCapture.isOpened())
  {
    std::cerr << "Unable to open the video : " << videoPath ;
    throw std::invalid_argument("Unable to open the video : "+videoPath);
  }
  // Grab frame 0, so we can call readImage.
  goToFrame(0);

  // load the calibration path
  _withIntrinsics = !calibPath.empty();
  if(_withIntrinsics)
    readCalibrationFromFile(calibPath, _camIntrinsics);
  
  _isInit = true;
}



bool VideoFeed::FeederImpl::readImage(image::Image<unsigned char> &imageGray,
                   cameras::Pinhole_Intrinsic_Radial_K3 &camIntrinsics,
                   std::string &mediaPath,
                   bool &hasIntrinsics)
{
  cv::Mat frame;
  const bool grabStatus = _videoCapture.retrieve(frame);

  if(!grabStatus || !frame.data)
  {
    return false;
  }
  
  if(frame.channels() == 3)
  {
    // convert to gray
    cv::Mat grey;
    cv::cvtColor(frame, grey, CV_BGR2GRAY);
    imageGray.resize(grey.cols, grey.rows);
    cv::cv2eigen(grey, imageGray);
//      std::cout << grey.channels() << " " << grey.rows << " " << grey.cols << std::endl;
//      std::cout << imageGray.Depth() << " " << imageGray.Height() << " " << imageGray.Width() << std::endl;
  }
  else
  {
    cv::cv2eigen(frame, imageGray);
  }

  hasIntrinsics = _withIntrinsics;
  if(_withIntrinsics)
    camIntrinsics = _camIntrinsics;

  mediaPath = _videoPath;
  return true;
}

std::size_t VideoFeed::FeederImpl::nbFrames() const
{
  if (!_videoCapture.isOpened())
    return 0;
  return _videoCapture.get(cv::CAP_PROP_FRAME_COUNT);
}

bool VideoFeed::FeederImpl::goToFrame(const unsigned int frame)
{
  if (!_videoCapture.isOpened())
  {
    std::cerr << "We cannot open the video file." << std::endl;
    return false;
  }
  
  if(frame > 0)
  {
    _videoCapture.set(cv::CAP_PROP_POS_FRAMES, frame);
    _videoCapture.grab();
    return true;
  }
  else
  {
    _videoCapture.set(cv::CAP_PROP_POS_FRAMES, 0);
    _videoCapture.grab();
    return false;
  }
}

bool VideoFeed::FeederImpl::goToNextFrame()
{
  return _videoCapture.grab();
}

/*******************************************************************************/
/*                                 VideoFeed                                   */
/*******************************************************************************/

VideoFeed::VideoFeed() : _feeder(new FeederImpl()) { }

VideoFeed::VideoFeed(const std::string &videoPath, const std::string &calibPath) 
  : _feeder(new FeederImpl(videoPath, calibPath))
{ }

bool VideoFeed::readImage(image::Image<unsigned char> &imageGray,
                     cameras::Pinhole_Intrinsic_Radial_K3 &camIntrinsics,
                     std::string &mediaPath,
                     bool &hasIntrinsics)
{
  return(_feeder->readImage(imageGray, camIntrinsics, mediaPath, hasIntrinsics));
}

std::size_t VideoFeed::nbFrames() const
{
  return _feeder->nbFrames();
}

bool VideoFeed::goToFrame(const unsigned int frame)
{
  return _feeder->goToFrame(frame);
}

bool VideoFeed::goToNextFrame()
{
  return _feeder->goToNextFrame();
}

bool VideoFeed::isInit() const {return(_feeder->isInit()); }

VideoFeed::~VideoFeed() { }

}//namespace dataio 
}//namespace openMVG


#endif //#if HAVE_OPENCV