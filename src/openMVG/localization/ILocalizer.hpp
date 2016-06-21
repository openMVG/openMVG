/* 
 * File:   ILocalizer.hpp
 * Author: sgaspari
 *
 * Created on November 18, 2015, 12:01 PM
 */

#pragma once

#include "LocalizationResult.hpp"

#include <openMVG/image/image_container.hpp>
#include <openMVG/features/image_describer.hpp>
#include <openMVG/cameras/Camera_Pinhole_Radial.hpp>

namespace openMVG {
namespace localization {

struct LocalizerParameters
{
  LocalizerParameters() :
  _visualDebug(""),
  _refineIntrinsics(false),
  _fDistRatio(0.8),
  _featurePreset(features::EDESCRIBER_PRESET::ULTRA_PRESET),
  _errorMax(std::numeric_limits<double>::max()) { }

  std::string _visualDebug;          //< enable visual debugging options
  bool _refineIntrinsics;     //< whether or not the Intrinsics of the query camera has to be refined
  float _fDistRatio;          //< the ratio distance to use when matching feature with the ratio test
  features::EDESCRIBER_PRESET _featurePreset; //< the preset to use for feature extraction of the query image
  double _errorMax;  
};

class ILocalizer
{
public:
    ILocalizer() : _isInit(false) { };
    
    bool isInit() {return _isInit;}
    
    const sfm::SfM_Data& getSfMData() const {return _sfm_data; }
    
    /**
   * @brief Localize one image
   * 
   * @param[in] imageGrey The input greyscale image.
   * @param[in] param The parameters for the localization.
   * @param[in] useInputIntrinsics Uses the \p queryIntrinsics as known calibration.
   * @param[in,out] queryIntrinsics Intrinsic parameters of the camera, they are used if the
   * flag useInputIntrinsics is set to true, otherwise they are estimated from the correspondences.
   * @param[out] localizationResult The localization result containing the pose and the associations.
   * @param[in] imagePath Optional complete path to the image, used only for debugging purposes.
   * @return  true if the image has been successfully localized.
   */
  virtual bool localize(const image::Image<unsigned char> & imageGrey,
                        const LocalizerParameters *param,
                        bool useInputIntrinsics,
                        cameras::Pinhole_Intrinsic_Radial_K3 &queryIntrinsics,
                        LocalizationResult & localizationResult,
                        const std::string& imagePath = std::string()) = 0;

  virtual bool localize(const std::unique_ptr<features::Regions> &queryRegions,
                        const std::pair<std::size_t, std::size_t> &imageSize,
                        const LocalizerParameters *param,
                        bool useInputIntrinsics,
                        cameras::Pinhole_Intrinsic_Radial_K3 &queryIntrinsics,
                        LocalizationResult & localizationResult,
                        const std::string& imagePath = std::string()) = 0;
    
  virtual bool localizeRig(const std::vector<image::Image<unsigned char> > & vec_imageGrey,
                           const LocalizerParameters *param,
                           std::vector<cameras::Pinhole_Intrinsic_Radial_K3 > &vec_queryIntrinsics,
                           const std::vector<geometry::Pose3 > &vec_subPoses,
                           geometry::Pose3 rigPose)=0;
   
  virtual ~ILocalizer( ) { } ;
protected:
  bool _isInit;
  sfm::SfM_Data _sfm_data;

};

} //namespace openMVG 
} //namespace localization 

