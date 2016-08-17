/* 
 * File:   FeedProvider.cpp
 * Author: sgaspari
 * 
 * Created on September 28, 2015, 6:48 PM
 */

#include "FeedProvider.hpp"
#include "ImageFeed.hpp"
#include "VideoFeed.hpp"

#include <boost/filesystem.hpp>

#include <exception>
#include <iostream>

namespace openMVG{
namespace dataio{

FeedProvider::FeedProvider(const std::string &feedPath, const std::string &calibPath) : _isVideo(false)
{
  namespace bf = boost::filesystem;
  if(feedPath.empty())
  {
    throw std::invalid_argument("Empty filepath.");
  }
  if(bf::is_regular_file(bf::path(feedPath))) 
  {
    // Image or video file
    const std::string extension = bf::path(feedPath).extension().string();
    if(ImageFeed::isSupported(extension))
    {
      _feeder.reset(new ImageFeed(feedPath, calibPath));
    }
    else 
    {
#if HAVE_OPENCV
      // let's try it with a video
      _feeder.reset(new VideoFeed(feedPath, calibPath));
      _isVideo = true;
#else
      throw std::invalid_argument("Unsupported mode! If you intended to use a video"
                                  " please add OpenCV support");
#endif
    }
  }
  // parent_path() returns "/foo/bar/" when input path equals to "/foo/bar/"
  // if the user just gives the relative path as "bar", throws invalid argument exception.
  else if(bf::is_directory(bf::path(feedPath)) || bf::is_directory(bf::path(feedPath).parent_path()))
  {
    // Folder or sequence of images
    _feeder.reset(new ImageFeed(feedPath, calibPath));
  }
  else
  {
    throw std::invalid_argument(std::string("Input filepath not supported: ") + feedPath);
  }
}

bool FeedProvider::readImage(image::Image<unsigned char> &imageGray,
      cameras::Pinhole_Intrinsic_Radial_K3 &camIntrinsics,
      std::string &mediaPath,
      bool &hasIntrinsics)
{
  return(_feeder->readImage(imageGray, camIntrinsics, mediaPath, hasIntrinsics));
}
  
std::size_t FeedProvider::nbFrames() const
{
  return _feeder->nbFrames();
}

bool FeedProvider::goToFrame(const unsigned int frame)
{
  return _feeder->goToFrame(frame);
}

bool FeedProvider::goToNextFrame()
{
  return _feeder->goToNextFrame();
}

bool FeedProvider::isInit() const
{
  return(_feeder->isInit());
}

FeedProvider::~FeedProvider( ) { }

}//namespace dataio 
}//namespace openMVG
