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

FeedProvider::FeedProvider(const std::string &feedPath, const std::string &calibPath)
{
  namespace bf = boost::filesystem;
  if(!bf::exists(feedPath))
  {
    throw std::invalid_argument(feedPath+" does not exist!");
  }
  
  if(bf::is_regular_file(bf::path(feedPath))) 
  {
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
#else
      throw std::invalid_argument("Unsupported mode! If you intended to use a video"
              " please add OpenCV support");
#endif
    }
  }
  else if(bf::is_directory(bf::path(feedPath)))
  {
    std::cout << "directory\n"; 
    _feeder.reset(new ImageFeed(feedPath, calibPath));
  }
  else
  {
    // throw something
    throw std::invalid_argument("File or mode not (yet) supported");
  }
}

  
bool FeedProvider::next(image::Image<unsigned char> &imageGray,
      cameras::Pinhole_Intrinsic_Radial_K3 &camIntrinsics,
      std::string &mediaPath,
      bool &hasIntrinsics)
{
  return(_feeder->next(imageGray, camIntrinsics, mediaPath, hasIntrinsics));
}

bool FeedProvider::isInit() const
{
  return(_feeder->isInit());
}

FeedProvider::~FeedProvider( ) { }


}//namespace dataio 
}//namespace openMVG