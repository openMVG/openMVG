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
    void serialize(Archive & ar)
    {
      ar(options_, eAkazeDescriptor_);
    }

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
    switch(preset)
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
  void serialize(Archive & ar)
  {
    ar(
     cereal::make_nvp("params", params_),
     cereal::make_nvp("bOrientation", bOrientation_));
  }

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

  std::unique_ptr<Regions_type> Describe(
    const image::Image<unsigned char>& image,
    const image::Image<unsigned char>* mask = nullptr
  );

  std::unique_ptr<Regions_type> Allocate() const
  {
    return std::unique_ptr<Regions_type>(new Regions_type);
  }

private:
  std::unique_ptr<Regions> DescribeImpl(
    const image::Image<unsigned char>& image,
    const image::Image<unsigned char>* mask = nullptr
  ) override
  {
    return Describe(image, mask);
  }

  virtual std::unique_ptr<Regions> AllocateImpl() const override
  {
    return Allocate();
  }
};

class AKAZE_Image_describer_LIOP : public AKAZE_Image_describer
{

public:
  using Regions_type = AKAZE_Liop_Regions;

  AKAZE_Image_describer_LIOP(
    const Params& params = Params(),
    bool bOrientation = true
  ):AKAZE_Image_describer(params, bOrientation) { }

  std::unique_ptr<Regions_type> Describe(
    const image::Image<unsigned char>& image,
    const image::Image<unsigned char>* mask = nullptr
  );

  std::unique_ptr<Regions_type> Allocate() const
  {
    return std::unique_ptr<Regions_type>(new Regions_type);
  }

protected:
  std::unique_ptr<Regions> DescribeImpl(
    const image::Image<unsigned char>& image,
    const image::Image<unsigned char>* mask = nullptr
  ) override
  {
    return Describe(image, mask);
  }

  virtual std::unique_ptr<Regions> AllocateImpl() const override
  {
    return Allocate();
  }
};

class AKAZE_Image_describer_MLDB : public AKAZE_Image_describer
{
public:
  using Regions_type = AKAZE_Binary_Regions;

  AKAZE_Image_describer_MLDB(
    const Params& params = Params(),
    bool bOrientation = true
  ):AKAZE_Image_describer(params, bOrientation) { }

  std::unique_ptr<Regions_type> Describe(
    const image::Image<unsigned char>& image,
    const image::Image<unsigned char>* mask = nullptr
  );

  std::unique_ptr<Regions_type> Allocate() const
  {
    return std::unique_ptr<Regions_type>(new Regions_type);
  }

protected:
  float GetfDescFactor() const override
  {
    return 11.f*sqrtf(2.f);
  }

  std::unique_ptr<Regions> DescribeImpl(
    const image::Image<unsigned char>& image,
    const image::Image<unsigned char>* mask = nullptr
  ) override
  {
    return Describe(image, mask);
  }

  virtual std::unique_ptr<Regions> AllocateImpl() const override
  {
    return Allocate();
  }
};

} // namespace features
} // namespace openMVG

#include <cereal/types/polymorphic.hpp>
#include <cereal/archives/json.hpp>
CEREAL_REGISTER_TYPE_WITH_NAME(openMVG::features::AKAZE_Image_describer, "AKAZE_Image_describer");
CEREAL_REGISTER_POLYMORPHIC_RELATION(openMVG::features::Image_describer, openMVG::features::AKAZE_Image_describer)

#endif // OPENMVG_FEATURES_AKAZE_IMAGE_DESCRIBER_HPP
