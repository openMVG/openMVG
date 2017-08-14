// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// The <cereal/archives> headers are special and must be included first.
#include <cereal/archives/json.hpp>

#include "openMVG/features/image_describer_akaze_io.hpp"
#include "openMVG/features/sift/SIFT_Anatomy_Image_Describer_io.hpp"
#include "openMVG/features/regions_factory_io.hpp"

#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include "testing/testing.h"

#include <fstream>

using namespace openMVG;
using namespace openMVG::features;
using namespace std;

using std::string;

bool SaveAndLoad
(
  std::unique_ptr<Image_describer> & image_describer
)
{
  const std::string sImage_describer = stlplus::create_filespec("./", "image_describer", "json");
  {
    std::ofstream stream(sImage_describer.c_str());
    if (!stream.is_open())
      return false;

    try
    {
      cereal::JSONOutputArchive archive(stream);
      archive(cereal::make_nvp("image_describer", image_describer));
      auto regionsType = image_describer->Allocate();
      archive(cereal::make_nvp("regions_type", regionsType));
    }
    catch (const cereal::Exception & e)
    {
      std::cerr << e.what() << std::endl
        << "Cannot dynamically allocate the Image_describer interface." << std::endl;
      return false;
    }
  }

  image_describer.reset();

  {
    // Dynamically load the image_describer from the file (will restore old used settings)
    std::ifstream stream(sImage_describer.c_str());
    if (!stream.is_open())
      return false;
    try
    {
      cereal::JSONInputArchive archive(stream);
      archive(cereal::make_nvp("image_describer", image_describer));
    }
    catch (const cereal::Exception & e)
    {
      std::cerr << e.what() << std::endl
        << "Cannot dynamically allocate the Image_describer interface." << std::endl;
      return false;
    }
  }
  return true;
}

TEST(Image_describer_akaze_surf, IO)
{
  std::unique_ptr<Image_describer> image_describer = AKAZE_Image_describer::create
    (AKAZE_Image_describer::Params(AKAZE::Params(), AKAZE_MSURF));

  EXPECT_TRUE(SaveAndLoad(image_describer));
}

TEST(Image_describer_akaze_mldb, IO)
{
  std::unique_ptr<Image_describer> image_describer = AKAZE_Image_describer::create
    (AKAZE_Image_describer::Params(AKAZE::Params(), AKAZE_MLDB));

  EXPECT_TRUE(SaveAndLoad(image_describer));
}

TEST(Image_describer_sift_anatomy, IO)
{
  std::unique_ptr<Image_describer> image_describer(new SIFT_Anatomy_Image_describer(SIFT_Anatomy_Image_describer::Params()));

  EXPECT_TRUE(SaveAndLoad(image_describer));
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
