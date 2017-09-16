// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Romuald PERROT, Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/features/feature.hpp"
#include "openMVG/features/mser/mser.hpp"
#include "openMVG/features/mser/mser_region.hpp"
#include "openMVG/image/image_io.hpp"

#include "testing/testing.h"


using namespace openMVG;
using namespace image;
using namespace features;
using namespace MSER;

TEST( MSER , Extraction )
{
  Image<unsigned char> image , outimg;
  const std::string png_filename = std::string(THIS_SOURCE_DIR) + "/image_test/lena.png";
  EXPECT_TRUE(ReadImage(png_filename.c_str(), &image));

  // Inverted image
  Image<unsigned char> image4( 255 - image.array() );

  std::vector< MSERRegion > regs;

  MSERExtractor extr4( 2 , 0.0005, 0.1 , 0.5 , 0.5 , MSERExtractor::MSER_4_CONNECTIVITY );
  extr4.Extract( image4 , regs );
  MSERExtractor extr8( 2 , 0.0005 , 0.1 , 0.5 , 0.5 , MSERExtractor::MSER_8_CONNECTIVITY );
  extr8.Extract( image , regs );

  // Export ellipse
  outimg = image;
  for (size_t i = 0; i < regs.size(); ++i )
  {
    double ex , ey;
    double Mx, My;
    double mx , my;
    double Ml , ml;
    regs[i].FitEllipse( ex , ey , Mx , My , mx , my , Ml , ml );

    for (double t = 0; t < 2.0 * 3.141592; t += 0.001 )
    {
      const int x = ex + ( cos( t ) * Mx * Ml + sin( t ) * mx * ml ) * 2.0 + 0.5;
      const int y = ey + ( cos( t ) * My * Ml + sin( t ) * my * ml ) * 2.0 + 0.5;

      if (x >= 0 && y >= 0 && x < image.Width() && y < image.Height() )
      {
        outimg( y , x ) = 255;
      }
    }

    /*
    double a, b, c;
    regs[i].FitEllipse( a, b, c );
    double x, y;
    regs[i].FitEllipse( x, y );

    const AffinePointFeature fp(x, y, a, b, c);
    DrawEllipse(fp.x(), fp.y(), fp.l1(), fp.l2(), 255, &outimg, fp.orientation());
    */
  }

  WriteImage( "outimg.png" , outimg );
}

TEST( MSER , Extraction_Synthetic )
{
  // Dark MSER detection
  {
    Image<unsigned char> image(40, 40, true, 254);
    image.block<10,10>(15,15).setConstant(1);
    WriteImage( "in_1.png" , image );

    std::vector< MSERRegion > regs;
    // Dark regions
    MSERExtractor extr( 2 , 0.0005 , 0.1 , 0.5 , 0.5 , MSERExtractor::MSER_4_CONNECTIVITY );
    extr.Extract( image , regs );

    EXPECT_EQ(1, regs.size());
  }

  // Bright MSER detection
  {
    Image<unsigned char> image(40, 40, true, 1);
    image.block<10,10>(15,15).setConstant(127);
    WriteImage( "in_2.png" , image );

    // Dark regions
    std::vector< MSERRegion > regs;
    MSERExtractor extr( 2 , 0.0005 , 0.1 , 0.5 , 0.5 , MSERExtractor::MSER_4_CONNECTIVITY );
    extr.Extract( Image<unsigned char>( 255 - image.array() ) , regs );

    EXPECT_EQ(1, regs.size());
  }
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
