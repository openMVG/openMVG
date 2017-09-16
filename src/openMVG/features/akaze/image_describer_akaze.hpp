// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_FEATURES_AKAZE_IMAGE_DESCRIBER_HPP
#define OPENMVG_FEATURES_AKAZE_IMAGE_DESCRIBER_HPP


#include "openMVG/features/akaze/AKAZE.hpp"
#include "openMVG/features/akaze/mldb_descriptor.hpp"
#include "openMVG/features/akaze/msurf_descriptor.hpp"
#include "openMVG/features/image_describer.hpp"
#include "openMVG/features/liop/liop_descriptor.hpp"
#include "openMVG/features/regions_factory.hpp"

#include <cereal/cereal.hpp>

#include <iostream>
#include <numeric>

namespace openMVG {
namespace features {


enum EAKAZE_DESCRIPTOR
{
  AKAZE_MSURF,
  AKAZE_LIOP,
  AKAZE_MLDB
};

class AKAZE_Image_describer : public Image_describer
{
public:

  struct Params
  {
    Params(
      const features::AKAZE::Params config = features::AKAZE::Params(),
      EAKAZE_DESCRIPTOR eAkazeDescriptor = AKAZE_MSURF
    ):options_(config), eAkazeDescriptor_(eAkazeDescriptor){}

    template<class Archive>
    void serialize(Archive & ar);

    // Parameters
    features::AKAZE::Params options_;
    EAKAZE_DESCRIPTOR eAkazeDescriptor_;
  };

  AKAZE_Image_describer(
    const Params & params = Params(),
    bool bOrientation = true
  ):Image_describer(), params_(params), bOrientation_(bOrientation) {}

  static std::unique_ptr<AKAZE_Image_describer> create(const Params& params, bool orientation = true);

  bool Set_configuration_preset(EDESCRIBER_PRESET preset) override
  {
    switch (preset)
    {
    case NORMAL_PRESET:
      params_.options_.fThreshold = features::AKAZE::Params().fThreshold;
    break;
    case HIGH_PRESET:
      params_.options_.fThreshold = features::AKAZE::Params().fThreshold/10.f;
    break;
    case ULTRA_PRESET:
      params_.options_.fThreshold = features::AKAZE::Params().fThreshold/100.f;
    break;
    default:
      return false;
    }
    return true;
  }

  template<class Archive>
  void serialize(Archive & ar);

protected:
  virtual float GetfDescFactor() const
  {
    return 10.f*sqrtf(2.f);
  }

  Params params_;
  bool bOrientation_;
};

class AKAZE_Image_describer_SURF : public AKAZE_Image_describer
{
public:

  using Regions_type = AKAZE_Float_Regions;

  AKAZE_Image_describer_SURF(
    const Params& params = Params(),
    bool bOrientation = true
  )
    :AKAZE_Image_describer(params, bOrientation) { }

  std::unique_ptr<Regions> Describe(
      const image::Image<unsigned char>& image,
      const image::Image<unsigned char>* mask = nullptr
  ) override
  {
    return Describe_AKAZE_SURF(image, mask);
  }

  std::unique_ptr<Regions> Allocate() const override
  {
    return std::unique_ptr<Regions_type>(new Regions_type);
  }

  std::unique_ptr<Regions_type> Describe_AKAZE_SURF(
    const image::Image<unsigned char>& image,
    const image::Image<unsigned char>* mask = nullptr
  );
};

class AKAZE_Image_describer_LIOP : public AKAZE_Image_describer {
public:
  using Regions_type = AKAZE_Liop_Regions;

  AKAZE_Image_describer_LIOP(
    const Params& params = Params(),
    bool bOrientation = true
  ):AKAZE_Image_describer(params, bOrientation) { }

  std::unique_ptr<Regions> Describe(
      const image::Image<unsigned char>& image,
      const image::Image<unsigned char>* mask = nullptr
  ) override
  {
    return Describe_AKAZE_LIOP(image, mask);
  }

  std::unique_ptr<Regions> Allocate() const override
  {
    return std::unique_ptr<Regions_type>(new Regions_type);
  }

  std::unique_ptr<Regions_type> Describe_AKAZE_LIOP(
    const image::Image<unsigned char>& image,
    const image::Image<unsigned char>* mask = nullptr
  );
};

class AKAZE_Image_describer_MLDB : public AKAZE_Image_describer
{
public:
  using Regions_type = AKAZE_Binary_Regions;

  AKAZE_Image_describer_MLDB(
    const Params& params = Params(),
    bool bOrientation = true
  ):AKAZE_Image_describer(params, bOrientation) { }

  std::unique_ptr<Regions> Describe(
    const image::Image<unsigned char>& image,
    const image::Image<unsigned char>* mask = nullptr
  ) override
  {
    return Describe_AKAZE_MLDB(image, mask);
  }

  std::unique_ptr<Regions> Allocate() const override
  {
    return std::unique_ptr<Regions_type>(new Regions_type);
  }

  std::unique_ptr<Regions_type> Describe_AKAZE_MLDB(
    const image::Image<unsigned char>& image,
    const image::Image<unsigned char>* mask = nullptr
  );

protected:
  float GetfDescFactor() const override
  {
    return 11.f*sqrtf(2.f);
  }
};

} // namespace features
} // namespace openMVG

#endif // OPENMVG_FEATURES_AKAZE_IMAGE_DESCRIBER_HPP
