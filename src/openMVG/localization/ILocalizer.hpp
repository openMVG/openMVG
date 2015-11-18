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
  _refineIntrinsics(false),
  _fDistRatio(0.6),
  _featurePreset(features::EDESCRIBER_PRESET::ULTRA_PRESET),
  _errorMax(std::numeric_limits<double>::max()) { }

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
    
    /**
    * @brief Just a wrapper around the different localization algorithm, the algorith
    * used to localized is chosen using \p param._algorithm
    * 
    * @param[in] imageGray The input greyscale image
    * @param[in] param The parameters for the localization
    * @param[in] useInputIntrinsics Uses the \p queryIntrinsics as known calibration
    * @param[in,out] queryIntrinsics Intrinsic parameters of the camera, they are used if the
    * flag useInputIntrinsics is set to true, otherwise they are estimated from the correspondences.
    * @param[out] pose The camera pose
    * @param[out] resection_data the 2D-3D correspondences used to compute the pose
    * @return true if the localization is successful
    */
    virtual bool localize(const image::Image<unsigned char> & imageGrey,
                 const LocalizerParameters &param,
                 bool useInputIntrinsics,
                 cameras::Pinhole_Intrinsic_Radial_K3 &queryIntrinsics,
                 LocalizationResult & localizationResult)=0;
   
    virtual ~ILocalizer( );
private:
  bool _isInit;

};

} //namespace openMVG 
} //namespace localization 

