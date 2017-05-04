// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_PATENTED_SIFT_SIFT_DESCRIBER_HPP
#define OPENMVG_PATENTED_SIFT_SIFT_DESCRIBER_HPP

#include "openMVG/features/image_describer.hpp"

namespace openMVG {
namespace features {

// Bibliography:
// [1] R. Arandjelović, A. Zisserman.
// Three things everyone should know to improve object retrieval. CVPR2012.

class SIFT_Image_describer : public Image_describer
{
public:

  struct Params
  {
    Params(
      int first_octave = 0,
      int num_octaves = 6,
      int num_scales = 3,
      float edge_threshold = 10.0f,
      float peak_threshold = 0.04f,
      bool root_sift = true
    );

    template<class Archive>
    void serialize( Archive & ar );

    // Parameters
    int _first_octave;      // Use original image, or perform an upscale if == -1
    int _num_octaves;       // Max octaves count
    int _num_scales;        // Scales per octave
    float _edge_threshold;  // Max ratio of Hessian eigenvalues
    float _peak_threshold;  // Min contrast
    bool _root_sift;        // see [1]
  };

  //--
  // Constructor
  //--
  SIFT_Image_describer
  (
    const Params params = Params(),
    bool bOrientation = true
  );

  ~SIFT_Image_describer();

  bool Set_configuration_preset(EDESCRIBER_PRESET preset) override;

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
  void serialize( Archive & ar );

private:
  Params _params;
  bool _bOrientation;
};

} // namespace features
} // namespace openMVG

#endif // OPENMVG_PATENTED_SIFT_SIFT_DESCRIBER_HPP
