/* 
 * File:   ImageFeed.hpp
 * Author: sgaspari
 *
 * Created on September 28, 2015, 10:31 AM
 */

#pragma once

#include "IFeed.hpp"

#include <string>
#include <memory>

namespace openMVG{
namespace dataio{

class ImageFeed : public IFeed
{
public:
  /**
   * @brief Empty constructor
   */	
  ImageFeed();
  
  /**
   * @brief Set up an image based feed from a choice of different sources:
   * 1) a directory containing images
   * 2) a json file containing a sfm data reconstruction (in that case \p calibPath is ignored)
   * 3) a txt file containing a list of images to use
   * 4) a regex for an image sequence
   * 
   * @param[in] imagePath The source of images, it could be a file (json, txt) or a directory.
   * @param[in] calibPath The source for the camera intrinsics common to each image. 
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
  ImageFeed(const std::string& imagePath, const std::string& calibPath);

  /**
   * @brief Provide a new image from the feed
   * 
   * @param[out] imageGray The new image from the feed.
   * @param[out] camIntrinsics The associated camera intrinsics.
   * @param[out] mediaPath The path to the current image.
   * @param[out] hasIntrinsics True if \p camIntrinsics is valid, otherwise there
   * is no intrinsics associated to \p imageGray.
   * @return True if there is a new image, false otherwise.
   */
  bool readImage(image::Image<unsigned char> &imageGray,
            cameras::Pinhole_Intrinsic_Radial_K3 &camIntrinsics,
            std::string &mediaPath,
            bool &hasIntrinsics);
  
  std::size_t nbFrames() const;
  
  bool goToFrame(const unsigned int frame);

  bool goToNextFrame();
  
  /**
   * @brief Return true if the feed is correctly initialized.
   * 
   * @return True if the feed is correctly initialized.
   */
  bool isInit() const;

  virtual ~ImageFeed( );
  
  /**
   * @brief For a given extension, return true if that file can be used as input
   * for the feed. ImageFeed supports .json, .txt, and the most common image files. 
   * 
   * @param extension The file extension to check in ".ext" format (case insensitive)
   * @return True if the file is supported.
   */
  static bool isSupported(const std::string &extension);
  
private:
  class FeederImpl;
  std::unique_ptr<FeederImpl> _imageFeed;
};

}//namespace dataio 
}//namespace openMVG


