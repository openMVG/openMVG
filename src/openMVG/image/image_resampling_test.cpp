// Copyright (c) 2015 Romuald Perrot.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/image/image.hpp"
#include "testing/testing.h"

#include <string>
#include <sstream>

using namespace openMVG;
using namespace openMVG::image;

TEST(Ressampling,SampleSamePosition)
{
  Image<unsigned char> image;
  const std::string png_filename = std::string(THIS_SOURCE_DIR) + "/image_test/lena.png";
  std::cout << png_filename << std::endl ;
  EXPECT_TRUE(ReadImage(png_filename.c_str(), &image));


  // Build sampling grid
  std::vector< std::pair< float , float > > sampling_grid ;
  sampling_grid.reserve(image.Width()*image.Height()) ;
  for( int i = 0 ; i < image.Height() ; ++i )
  {
    for( int j = 0 ; j < image.Width() ; ++j )
    {
      sampling_grid.emplace_back( i , j ) ;
    }
  }

  // Ressample image
  Sampler2d< SamplerLinear > sampler ;

  Image<unsigned char> imageOut ;

  GenericRessample( image , sampling_grid , image.Width() , image.Height() , sampler , imageOut ) ;

  const std::string out_filename = ("test_ressample_same.png");
  EXPECT_TRUE( WriteImage( out_filename.c_str(), imageOut) ) ;
}

// Iterative image rotations
// Allow to check if the sampling function have some signal loss.
template<typename SamplerT, typename ImageT>
bool ImageRotation(
  const ImageT & imageIn,
  const SamplerT & sampler,
  const std::string & samplerString)
{
  bool bOk = true;
  ImageT image = imageIn;

  const int nb_rot = 6 ;
  const float delta = ( 2.0 * 3.141592 ) / nb_rot ;

  const float middle_x = image.Width() / 2.0 ;
  const float middle_y = image.Height() / 2.0 ;

  // Rotate image then set starting image as source
  for( int id_rot = 0 ; id_rot < nb_rot ; ++id_rot )
  {
    // angle of rotation (negative because it's inverse transformation)
    const float cur_angle = delta ;

    const float cs = cosf( -cur_angle ) ;
    const float ss = sinf( -cur_angle ) ;

    std::vector< std::pair<float,float> > sampling_grid;
    sampling_grid.reserve(image.Width()*image.Height());
    // Compute sampling grid
    for( int i = 0 ; i < image.Height() ; ++i )
    {
      for( int j = 0 ; j < image.Width() ; ++j )
      {
        // Compute rotation of pixel (i,j) around center of image

        // Center pixel
        const float dx = static_cast<float>(j) - middle_x ;
        const float dy = static_cast<float>(i) - middle_y ;

        const float rotated_x = cs * dx - ss * dy ;
        const float rotated_y = ss * dx + cs * dy ;

        // Get back to original center
        const float cur_x = rotated_x + middle_x ;
        const float cur_y = rotated_y + middle_y ;

        sampling_grid.emplace_back( cur_y , cur_x) ;
      }
    }

    // Sample input image
    ImageT imageOut ;

    GenericRessample( image , sampling_grid , image.Width() , image.Height() , sampler , imageOut ) ;

    std::stringstream str ;
    str << "test_ressample_"<< samplerString <<"_rotate_" << id_rot << ".png" ;
    bOk &= WriteImage( str.str().c_str(), imageOut);
    image = imageOut ;
  }
  return bOk;
}

TEST(Ressampling,SampleRotate)
{
  Image<RGBColor> image;

  const std::string png_filename = std::string(THIS_SOURCE_DIR) + "/image_test/lena.png";
  std::cout << png_filename << std::endl ;
  EXPECT_TRUE(ReadImage(png_filename.c_str(), &image));

  EXPECT_TRUE(ImageRotation(image, Sampler2d< SamplerNearest >(), "SamplerNearest"));
  EXPECT_TRUE(ImageRotation(image, Sampler2d< SamplerLinear >(), "SamplerLinear"));
  EXPECT_TRUE(ImageRotation(image, Sampler2d< SamplerCubic >(), "SamplerCubic"));
  EXPECT_TRUE(ImageRotation(image, Sampler2d< SamplerSpline16 >(), "SamplerSpline16"));
  EXPECT_TRUE(ImageRotation(image, Sampler2d< SamplerSpline64 >(), "SamplerSpline64"));
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
