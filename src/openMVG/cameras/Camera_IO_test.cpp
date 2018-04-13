// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// The <cereal/archives> headers are special and must be included first.
#include <cereal/archives/json.hpp>

#include "openMVG/cameras/cameras_io.hpp"

#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <fstream>

using namespace openMVG;
using namespace openMVG::cameras;

#include "testing/testing.h"

using std::string;

TEST(Camera_IO_ceral, SaveRead) {

  const std::vector<EINTRINSIC> vec_camera_model_type =
    {
      PINHOLE_CAMERA,
      PINHOLE_CAMERA_RADIAL1, PINHOLE_CAMERA_RADIAL3,
      PINHOLE_CAMERA_BROWN,
      PINHOLE_CAMERA_FISHEYE
    };

  for (const auto cam_type : vec_camera_model_type)
  {
    std::shared_ptr<IntrinsicBase> intrinsic(nullptr);

    const int width = 200;
    const int height = 200;
    const double ppx = width / 2.0;
    const double ppy = height / 2.0;
    const double focal = 200;

    // Create the desired camera type
    switch (cam_type)
    {
    case PINHOLE_CAMERA:
      intrinsic = std::make_shared<Pinhole_Intrinsic>
        (width, height, focal, ppx, ppy);
      break;
    case PINHOLE_CAMERA_RADIAL1:
      intrinsic = std::make_shared<Pinhole_Intrinsic_Radial_K1>
        (width, height, focal, ppx, ppy);
      break;
    case PINHOLE_CAMERA_RADIAL3:
      intrinsic = std::make_shared<Pinhole_Intrinsic_Radial_K3>
        (width, height, focal, ppx, ppy);
      break;
    case PINHOLE_CAMERA_BROWN:
      intrinsic = std::make_shared<Pinhole_Intrinsic_Brown_T2>
        (width, height, focal, ppx, ppy);
      break;
    case PINHOLE_CAMERA_FISHEYE:
      intrinsic = std::make_shared<Pinhole_Intrinsic_Fisheye>
        (width, height, focal, ppx, ppy);
      break;
    }

    const std::string filename("camera_io.json");

    // Writing
    {
      std::ofstream stream(filename, std::ios::binary | std::ios::out);
      CHECK(stream.is_open());

      cereal::JSONOutputArchive archive(stream);
      archive(cereal::make_nvp("intrinsics", intrinsic));
    }
    // Reading
    {
      std::ifstream stream(filename, std::ios::binary | std::ios::in);
      CHECK(stream.is_open());

      cereal::JSONInputArchive archive(stream);
      archive(cereal::make_nvp("intrinsics", intrinsic));
    }
    CHECK(stlplus::file_delete(filename));
  }

}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
