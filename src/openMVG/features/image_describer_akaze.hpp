// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_FEATURES_AKAZE_IMAGE_DESCRIBER_HPP
#define OPENMVG_FEATURES_AKAZE_IMAGE_DESCRIBER_HPP

#include "openMVG/features/image_describer.hpp"
#include "openMVG/features/akaze/AKAZE.hpp"

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
    );

    template<class Archive>
    void serialize(Archive & ar);

    // Parameters
    features::AKAZE::Params options_;
    EAKAZE_DESCRIPTOR eAkazeDescriptor_;
  };

  AKAZE_Image_describer(
    const Params & params = std::move(Params()),
    bool bOrientation = true
  );


  bool Set_configuration_preset(EDESCRIBER_PRESET preset) override;;

  /**
  @brief Detect regions on the image and compute their attributes (description)
  @param image Image.
  @param regions The detected regions and attributes (the caller must delete the allocated data)
  @param mask 8-bit gray image for keypoint filtering (optional).
     Non-zero values depict the region of interest.
  */
  bool Describe
  (
    const image::Image<unsigned char>& image,
    std::unique_ptr<Regions> &regions,
    const image::Image<unsigned char> * mask = nullptr
  ) override;

  /// Allocate Regions type depending of the Image_describer
  void Allocate(std::unique_ptr<Regions> &regions) const override;

  template<class Archive>
  void serialize(Archive & ar);

private:
  Params params_;
  bool bOrientation_;
};

} // namespace features
} // namespace openMVG

#endif // OPENMVG_FEATURES_AKAZE_IMAGE_DESCRIBER_HPP
