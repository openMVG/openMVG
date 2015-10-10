/* 
 * File:   ImageFeed.cpp
 * Author: sgaspari
 * 
 * Created on September 28, 2015, 10:31 AM
 */

#include "ImageFeed.hpp"
#include <openMVG/sfm/sfm_data.hpp>
#include <openMVG/sfm/sfm_data_io.hpp>
#include <openMVG/image/image_io.hpp>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string/case_conv.hpp> 

#include <queue>
#include <iostream>
#include <fstream>
#include <exception>

namespace openMVG{
namespace dataio{

class ImageFeed::FeederImpl
{
public:
  
  static bool isSupported(const std::string &ext);
  
  FeederImpl() : _isInit(false) {}
  
  FeederImpl(const std::string imagePath, const std::string calibPath);

  bool next(image::Image<unsigned char> &imageGray, 
                     cameras::Pinhole_Intrinsic_Radial_K3 &camIntrinsics,
                     std::string &imageName,
                     bool &hasIntrinsics);
  
  bool isInit() const {return _isInit;} 
  
private:
  
  bool feedWithJson(image::Image<unsigned char> &imageGray, 
                     cameras::Pinhole_Intrinsic_Radial_K3 &camIntrinsics,
                     std::string &imageName,
                     bool &hasIntrinsics);
  
private:
  static const std::vector<std::string> supportedExtensions;
  
private:
  bool _isInit;
  bool _withCalibration;
  std::queue<std::string> _images;
  cameras::Pinhole_Intrinsic_Radial_K3 _camIntrinsics;
  
  bool _jsonMode = false;
  sfm::SfM_Data _sfmdata;
  sfm::Views::const_iterator _viewIterator;
};

const std::vector<std::string> ImageFeed::FeederImpl::supportedExtensions = { ".jpg", ".jpeg", ".png", ".ppm" };

bool ImageFeed::FeederImpl::isSupported(const std::string &ext)
{
  const auto start = FeederImpl::supportedExtensions.begin();
  const auto end = FeederImpl::supportedExtensions.end();
  return(std::find(start, end, boost::to_lower_copy(ext)) != end);
}


ImageFeed::FeederImpl::FeederImpl(const std::string imagePath, const std::string calibPath) 
: _isInit(false), _withCalibration(false)
{
  namespace bf = boost::filesystem;
//    std::cout << imagePath << std::endl;
  // if it is a json, calibPath is neglected
  if(bf::is_regular_file(imagePath))
  {
    const std::string ext = bf::path(imagePath).extension().string();
    // if it is a sfmdata.json
    if(ext == ".json")
    {
      // load the json
      _isInit = sfm::Load(_sfmdata, imagePath, sfm::ESfM_Data(sfm::ESfM_Data::VIEWS | sfm::ESfM_Data::INTRINSICS));
      _viewIterator = _sfmdata.GetViews().begin();
      _jsonMode = true;
    }
    // if it is an image file
    else if(FeederImpl::isSupported(ext))
    {
      _images.push(imagePath);
      _withCalibration = !calibPath.empty();
      _jsonMode = false;
      _isInit = true;
    }
    // if it is an image file
    else if(ext == ".txt")
    {
      // we expect a simple txt file with a list of path to images relative to the 
      // location of the txt file itself
      std::fstream fs(imagePath, std::ios::in);
      std::string line;
      // parse each line of the text file
      while(getline(fs, line))
      {
        // compose the file name as the base path of the inputPath and
        // the filename just read
        const std::string filename = (bf::path(imagePath).parent_path() / line).string();
        _images.push(filename);
      }
      // Close file
      fs.close();
      _withCalibration = !calibPath.empty();
      _jsonMode = false;
      _isInit = true;
    }
    else
    {
      // no other file format are supported
      throw std::invalid_argument("File or mode not yet implemented");
    }
  }
  else if(bf::is_directory(imagePath))
  {
    std::cout << "directory feedImage\n";
    // if it is a directory, list all the images and add them to the list
    bf::directory_iterator iterator(imagePath);
    for(; iterator != bf::directory_iterator(); ++iterator)
    {
      // get the extension of the current file to check whether it is an image
      const std::string ext = iterator->path().extension().string();
      if(FeederImpl::isSupported(ext))
      {
        _images.push(iterator->path().string());
      }
    }
    _withCalibration = !calibPath.empty();
    _jsonMode = false;
    _isInit = true;
  }
  else
  {
    throw std::invalid_argument("File or mode not yet implemented");
  }
  
  // last thing: if _withCalibration is true it means that it is not a json and
  // a path to a calibration file has been passed
  // then load the calibration
  if(_withCalibration)
  {
    // load the calibration from calibPath
    readCalibrationFromFile(calibPath, _camIntrinsics);
  }
}



bool ImageFeed::FeederImpl::next(image::Image<unsigned char> &imageGray, 
                   cameras::Pinhole_Intrinsic_Radial_K3 &camIntrinsics,
                   std::string &imageName,
                   bool &hasIntrinsics)
{
  if(!_isInit)
  {
    std::cerr << "Image feed is not initialized " << std::endl;
    return false;
  }

  // dealing with json mode
  if(_jsonMode)
  {
    return(feedWithJson(imageGray, camIntrinsics, imageName, hasIntrinsics));
  }
  else
  {
    if(_images.empty())
    {
      return false;
    }
    
    if(_withCalibration)
    {
      // get the calibration
      camIntrinsics = _camIntrinsics;
      hasIntrinsics = true;
    }
    else
    {
      hasIntrinsics = false;
    }
    imageName = _images.front();
    std::cout << imageName << std::endl;
    if (!image::ReadImage(imageName.c_str(), &imageGray))
    {
      std::cerr << "Error while opening image " << imageName << std::endl;
      return false;
    }
    _images.pop();
    return true;
  }
  return true;
}


bool ImageFeed::FeederImpl::feedWithJson(image::Image<unsigned char> &imageGray, 
                   cameras::Pinhole_Intrinsic_Radial_K3 &camIntrinsics,
                   std::string &imageName,
                   bool &hasIntrinsics)
{
  // if there are no more images to process
  if(_viewIterator == _sfmdata.GetViews().end())
  {
    return false;
  }

  namespace bf = boost::filesystem;

  // get the image
  const std::string rootPath = _sfmdata.s_root_path;
  const sfm::View *view = _viewIterator->second.get();
  imageName = (bf::path(rootPath) / bf::path(view->s_Img_path)).string();
  if (!image::ReadImage(imageName.c_str(), &imageGray))
  {
    std::cerr << "Error while opening image " << imageName << std::endl;
    return false;
  }
  // get the associated Intrinsics
  if((view->id_intrinsic == UndefinedIndexT) or (!_sfmdata.GetIntrinsics().count(view->id_intrinsic)))
  {
    std::cout << "Image "<< imageName << " does not have associated intrinsics" << std::endl;
    hasIntrinsics = false;
  }
  else
  {
    const cameras::IntrinsicBase * cam = _sfmdata.GetIntrinsics().at(view->id_intrinsic).get();
    if(cam->getType() != cameras::EINTRINSIC::PINHOLE_CAMERA_RADIAL3)
    {
      std::cerr << "Only Pinhole_Intrinsic_Radial_K3 is supported" << std::endl;
      hasIntrinsics = false;
    }
    else
    {
      const cameras::Pinhole_Intrinsic_Radial_K3 * intrinsics = dynamic_cast<const cameras::Pinhole_Intrinsic_Radial_K3*>(cam) ;

      // simply copy values
      camIntrinsics = *intrinsics;
      hasIntrinsics = true;
    }
  }
  ++_viewIterator;
  return true;
}




/*******************************************************************************/
/*                     ImageFeed                                               */
/*******************************************************************************/

ImageFeed::ImageFeed() : _imageFeed(new FeederImpl()) { }

ImageFeed::ImageFeed(const std::string imagePath, const std::string calibPath)  
    : _imageFeed( new FeederImpl(imagePath, calibPath) ) { }

bool ImageFeed::next(image::Image<unsigned char> &imageGray, 
                     cameras::Pinhole_Intrinsic_Radial_K3 &camIntrinsics,
                     std::string &mediaPath,
                     bool &hasIntrinsics)
{
  return(_imageFeed->next(imageGray, camIntrinsics, mediaPath, hasIntrinsics));
}

bool ImageFeed::isInit() const
{
  return(_imageFeed->isInit());
}

bool ImageFeed::isSupported(const std::string &extension)
{
  std::string ext = boost::to_lower_copy(extension);
  if(ext == ".json" or ext == ".txt")
  {
    return true;
  }
  else
  {
    return FeederImpl::isSupported(ext);
  }
}

ImageFeed::~ImageFeed() { }

}//namespace dataio 
}//namespace openMVG