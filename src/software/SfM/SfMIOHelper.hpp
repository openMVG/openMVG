
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_IO_H
#define OPENMVG_SFM_IO_H

#include "openMVG/numeric/numeric.h"
#include "openMVG/split/split.hpp"

#include <fstream>
#include <iterator>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "JsonBox.h"

namespace openMVG{
namespace SfMIO{

struct CameraInfo
{
  std::string m_sImageName;
  size_t m_intrinsicId;
};

struct IntrinsicCameraInfo
{
  size_t m_w, m_h;
  float m_focal;
  Mat3 m_K;
  bool m_bKnownIntrinsic; // true if 11 or 6, else false
  std::string m_sCameraMaker, m_sCameraModel;

  bool operator() (IntrinsicCameraInfo const &ci1, IntrinsicCameraInfo const &ci2)const
  {
    bool bequal = false;
    if ( ci1.m_sCameraMaker.compare("") != 0  && ci1.m_sCameraModel.compare("") != 0 )
    {
      if ( ci1.m_sCameraMaker.compare(ci2.m_sCameraMaker) == 0
          && ci1.m_sCameraModel.compare(ci2.m_sCameraModel) == 0
          && ci1.m_w == ci2.m_w
          && ci1.m_h == ci2.m_h
          && ci1.m_focal == ci2.m_focal )
      {
        bequal = true;
      }
      else
      {
        if(m_bKnownIntrinsic)
          bequal = ci1.m_K == ci2.m_K;
      }
    }
    return !bequal;
  }
};



// Load an image file list
// One basename per line.
// It handle different scenario based on the intrinsic info of the tested image
// - a camera without exif data
// - a camera with exif data found in the database
// - a camera with exif data not found in the database
// - a camera with known intrinsic
static bool loadImageList( std::vector<CameraInfo> & vec_camImageName,
                           std::vector<IntrinsicCameraInfo> & vec_focalGroup,
                           std::string sFileName,
                           bool bVerbose = true )
{
  typedef std::set<IntrinsicCameraInfo, IntrinsicCameraInfo> setIntrinsicCameraInfo;
  setIntrinsicCameraInfo set_focalGroup;

  JsonBox::Value imageParams;
  imageParams.loadFromFile(sFileName);
  JsonBox::Array images;
  images = imageParams["images"].getArray();
  for(int i = 0; i < images.size(); i++)
  {
    IntrinsicCameraInfo intrinsicCamInfo;
    intrinsicCamInfo.m_w = images[i]["width"].getInt();
    intrinsicCamInfo.m_h = images[i]["height"].getInt();

    //a camera without exif
    intrinsicCamInfo.m_focal = -1;
    intrinsicCamInfo.m_bKnownIntrinsic = false;
    intrinsicCamInfo.m_sCameraMaker = "";
    intrinsicCamInfo.m_sCameraModel = "";

    //known camera make/model
    if(images[i]["camera"]["name"].isString())
      intrinsicCamInfo.m_sCameraMaker = images[i]["camera"]["name"].getString();
    if(images[i]["camera"]["model"].isString())
      intrinsicCamInfo.m_sCameraModel = images[i]["camera"]["model"].getString();

    //known focal
    if(images[i]["camera"]["focal"].getDouble() > 0.0)
    {
      intrinsicCamInfo.m_bKnownIntrinsic = true;
      intrinsicCamInfo.m_focal = images[i]["camera"]["focal"].getDouble();
      Mat3 K;
      K << intrinsicCamInfo.m_focal, 0, intrinsicCamInfo.m_w / 2,
          0, intrinsicCamInfo.m_focal, intrinsicCamInfo.m_h / 2,
          0, 0, 1;
      intrinsicCamInfo.m_K = K;
    }

    //if you knew K explicitly you could set it here too...
    //yes, storing the whole K is redundant since it looks like
    // f 0 w/2
    // 0 f h/2
    // 0 0 1
    if(images[i]["camera"]["K"].isArray())
    {
      Mat3 K;
      K << images[i]["camera"]["K"][size_t(0)].getDouble(),
           images[i]["camera"]["K"][size_t(1)].getDouble(),
           images[i]["camera"]["K"][size_t(2)].getDouble(),
           images[i]["camera"]["K"][size_t(3)].getDouble(),
           images[i]["camera"]["K"][size_t(4)].getDouble(),
           images[i]["camera"]["K"][size_t(5)].getDouble(),
           images[i]["camera"]["K"][size_t(6)].getDouble(),
           images[i]["camera"]["K"][size_t(7)].getDouble(),
           images[i]["camera"]["K"][size_t(8)].getDouble();
      intrinsicCamInfo.m_K = K;
    }

    // Setup Intrinsics...
    std::pair<setIntrinsicCameraInfo::iterator, bool> ret = set_focalGroup.insert(intrinsicCamInfo);
    if ( ret.second )
    {
      vec_focalGroup.push_back(intrinsicCamInfo);
    }
    size_t id = std::distance( ret.first, set_focalGroup.end()) - 1;
    CameraInfo camInfo;
    camInfo.m_sImageName = images[i]["filename"].getString();
    camInfo.m_intrinsicId = id;
    vec_camImageName.push_back(camInfo);
  }

  return !(vec_camImageName.empty());
}

} // namespace SfMIO
} // namespace openMVG

#endif // OPENMVG_SFM_INCREMENTAL_ENGINE_H

