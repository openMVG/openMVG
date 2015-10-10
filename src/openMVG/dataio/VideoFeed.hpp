/* 
 * File:   VideoFeed.hpp
 * Author: sgaspari
 *
 * Created on September 28, 2015, 10:35 AM
 */

#if HAVE_OPENCV

#pragma once

#include "IFeed.hpp"

#include <string>
#include <memory>

namespace openMVG{
namespace dataio{

class VideoFeed : public IFeed
{
public:
  VideoFeed();

  /**
   * @brief Set up an image feed from a video
   * 
   * @param[in] imagePath The video source.
   * @param[in] calibPath The source for the camera intrinsics. 
   * The format for the file is
   * int #image width
   * int #image height
   * double #focal
   * double #ppx principal point x-coord
   * double #ppy principal point y-coord
   * double #k0
   * double #k1
   * double #k2
   * 
   * @see readCalibrationFromFile()
   */  
  VideoFeed(const std::string &videoPath, const std::string &calibPath);

  /**
   * @brief Provide a new image from the feed
   * 
   * @param[out] imageGray The new image from the feed.
   * @param[out] camIntrinsics The associated camera intrinsics.
   * @param[out] mediaPath The original video path.
   * @param[out] hasIntrinsics True if \p camIntrinsics is valid, otherwise there
   * is no intrinsics associated to \p imageGray.
   * @return True if there is a new image, false otherwise.
   */
  bool next(image::Image<unsigned char> &imageGray,
            cameras::Pinhole_Intrinsic_Radial_K3 &camIntrinsics,
            std::string &mediaPath,
            bool &hasIntrinsics);
  /**
   * @brief Return true if the feed is correctly initialized.
   * 
   * @return True if the feed is correctly initialized.
   */  
  bool isInit() const;

  virtual ~VideoFeed( );
  
private:
  class FeederImpl;
  std::unique_ptr<FeederImpl> _feeder;
};

}//namespace dataio 
}//namespace openMVG

#endif //#if HAVE_OPENCV