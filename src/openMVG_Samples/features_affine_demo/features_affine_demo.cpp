// Copyright (c) 2015 Pierre MOULON, Romuald PERROT.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/image/image.hpp"
#include "openMVG/features/features.hpp"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "third_party/vectorGraphics/svgDrawer.hpp"
#include "third_party/cmdLine/cmdLine.h"

#include <string>
#include <iostream>

#include <unsupported/Eigen/MatrixFunctions>

using namespace openMVG;
using namespace openMVG::image;
using namespace openMVG::features;
using namespace svg;
using namespace std;

/**
 ** Normalize a patch to a specified size, given an ellipse
 ** Ie: build a patch
 **/
template <typename Image>
void NormalizePatch
(
  const Image & src_img ,
  const AffinePointFeature & feat ,
  const int patch_size ,
  Image & out_patch
)
{
  // Mapping function
  Eigen::Matrix<double,2,2> A;
  A << feat.a(), feat.b(),
       feat.b(), feat.c();

  // Inverse square root
  A = A.pow( -0.5 ) ;

  const float sc = 2.f * 3.f / static_cast<float>(patch_size);
  A = A * sc ;

  const float half_width = static_cast<float>( patch_size ) / 2.f ;

  // Compute sampling grid
  std::vector< std::pair<float,float> > sampling_grid ;
  sampling_grid.reserve(patch_size*patch_size);
  for( int i = 0 ; i < patch_size ; ++i )
  {
    for( int j = 0 ; j < patch_size ; ++j )
    {
      // Apply transformation relative to the center of the patch (assume origin at 0,0 then map to (x,y) )
      Vec2 pos;
      pos << static_cast<float>( j ) - half_width, static_cast<float>( i ) - half_width ;
      // Map (ie: ellipse transform)
      const Vec2 affineAdapted = A * pos;

      sampling_grid.emplace_back( affineAdapted(1) + feat.y() , affineAdapted(0) + feat.x() );
    }
  }

  Sampler2d< SamplerLinear > sampler ;

  // Sample input image to generate patch
  GenericRessample(
    src_img , sampling_grid ,
    patch_size , patch_size ,
    sampler ,
    out_patch ) ;
}

void Extract_MSER
(
  const Image<unsigned char> & img,
  std::vector<features::AffinePointFeature> & feats_dark,
  std::vector<features::AffinePointFeature> & feats_bright
)
{
  using namespace openMVG::features::MSER;

  //-- Extract Bight MSER
  {
    // Inverted image
    Image<unsigned char> image4( 255 - img.array() );
    std::vector< MSERRegion > regs;
    MSERExtractor extr4( 2 , 0.0005 , 0.1 , 0.5 , 0.5 , MSERExtractor::MSER_4_CONNECTIVITY ) ;
    extr4.Extract( image4 , regs ) ;
    for( int i = 0 ; i < regs.size() ; ++i )
    {
      double a, b, c;
      regs[i].FitEllipse( a, b, c );
      double x, y;
      regs[i].FitEllipse( x, y );
      feats_bright.emplace_back(x, y, a, b, c);
    }
  }

  //-- Extract Dark MSER
  {
    std::vector< MSERRegion > regs;
    MSERExtractor extr8( 2 , 0.0005 , 0.1 , 0.5 , 0.5 , MSERExtractor::MSER_8_CONNECTIVITY ) ;
    extr8.Extract( img , regs ) ;
    for( int i = 0 ; i < regs.size() ; ++i )
    {
      double a, b, c;
      regs[i].FitEllipse( a, b, c );
      double x, y;
      regs[i].FitEllipse( x, y );
      feats_dark.emplace_back(x, y, a, b, c);
    }
  }
}

void Extract_TBMR
(
  const Image<unsigned char> & img,
  std::vector<features::AffinePointFeature> & feats_dark,
  std::vector<features::AffinePointFeature> & feats_bright
)
{
  tbmr::Extract_tbmr (img, feats_bright, std::less<uint8_t> (), 30);
  tbmr::Extract_tbmr (img, feats_dark, std::greater<uint8_t> (), 30);
}

int main(int argc, char **argv)
{
  std::string sAffine_Detector_Method = "TBMR";

  CmdLine cmd;
  cmd.add(make_switch('P', "PATCH"));
  cmd.add(make_option('d', sAffine_Detector_Method, "detector") );

  std::cout
    << "TBMR Demo:\n"
    << " Show detected Affine regions as ellipses,\n"
    << " -[P] in the command line exports square normalized patches for each ellipses.\n"
    << " -[d|detector] TBMR|MSER Detect TBMR or MSER affine regions."
    << std::endl;

  try {
    cmd.process(argc, argv);
  } catch(const std::string& s) {
    std::cerr << s << std::endl;
    return EXIT_FAILURE;
  }
  const std::string sInputDir =
    stlplus::folder_up(string(THIS_SOURCE_DIR)) + "/imageData/SceauxCastle/";
  const std::string jpg_filename = sInputDir + "100_7101.jpg";

  Image<unsigned char> image;
  ReadImage(jpg_filename.c_str(), &image);

  std::vector<features::AffinePointFeature> feats_dark, feats_bright;
  if (sAffine_Detector_Method == "MSER")
  {
    Extract_MSER(image, feats_dark, feats_bright);
  }
  else if (sAffine_Detector_Method == "TBMR")
  {
    Extract_TBMR(image, feats_dark, feats_bright);
  }
  else
  {
    std::cerr << "Invalid Affine detector type." << std::endl;
    return EXIT_FAILURE;
  }

  //-- Affine Detector demo:
  {
    std::cout << "#detected BRIGHT " << sAffine_Detector_Method << ": " << feats_bright.size() << std::endl;

    // Display extracted Region ellipses:
    Image<unsigned char> Icpy (image);
    for( size_t i = 0; i < feats_bright.size(); ++i)
    {
      const AffinePointFeature & fp = feats_bright[i];
      DrawEllipse(fp.x(), fp.y(), fp.l1(), fp.l2(), 255, &Icpy, fp.orientation());
      if (cmd.used('P'))
      {
        //-- Ellipse to square 41x41 patch normalization
        Image<unsigned char> patch;
        NormalizePatch( Icpy , fp , 41 , patch );
        std::stringstream str;
        str << "Patch_" << i << ".png";
        WriteImage( str.str().c_str() , patch );
      }
    }
    std::ostringstream os;
    os << sAffine_Detector_Method << "_BRIGHT_features.jpg";
    WriteImage(os.str().c_str(), Icpy);

    std::cout << "#detected DARK " << sAffine_Detector_Method << ": " << feats_dark.size() << std::endl;

    // Display extracted Region ellipses:
    Icpy = image;
    for( size_t i = 0; i < feats_dark.size(); ++i)
    {
      const AffinePointFeature & fp = feats_dark[i];
      DrawEllipse(fp.x(), fp.y(), fp.l1(), fp.l2(), 255, &Icpy, fp.orientation());
    }
    os.str("");
    os << sAffine_Detector_Method << "_DARK_features.jpg";
    WriteImage(os.str().c_str(), Icpy);
  }

  return EXIT_SUCCESS;
}

