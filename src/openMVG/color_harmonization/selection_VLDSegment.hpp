// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2014 openMVG authors.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_COLOR_HARMONIZATION_SELECTION_VLDSEGMENT_HPP
#define OPENMVG_COLOR_HARMONIZATION_SELECTION_VLDSEGMENT_HPP

#include <string>
#include <vector>

#include "openMVG/color_harmonization/selection_interface.hpp"
#include "openMVG/matching/kvld/kvld.h"
#include "openMVG/matching/kvld/kvld_draw.h"

namespace openMVG {
namespace color_harmonization {

class commonDataByPair_VLDSegment  : public commonDataByPair
{
  public:
  commonDataByPair_VLDSegment( const std::string & sLeftImage,
                               const std::string & sRightImage,
                               const std::vector<matching::IndMatch>& vec_PutativeMatches,
                               const std::vector<features::SIOPointFeature >& vec_featsL,
                               const std::vector<features::SIOPointFeature >& vec_featsR):
           commonDataByPair( sLeftImage, sRightImage ),
           _vec_featsL( vec_featsL ), _vec_featsR( vec_featsR ),
           _vec_PutativeMatches( vec_PutativeMatches )
  {}

  ~commonDataByPair_VLDSegment() override = default;

  /**
   * Put masks to white, images are conserved
   *
   * \param[out] maskLeft Mask of the left image (initialized to corresponding image size).
   * \param[out] maskRight  Mask of the right image (initialized to corresponding image size).
   *
   * \return True.
   */
  bool computeMask(
    image::Image<unsigned char> & maskLeft,
    image::Image<unsigned char> & maskRight ) override
  {
    std::vector<matching::IndMatch> vec_KVLDMatches;

    image::Image<unsigned char> imageL, imageR;
    image::ReadImage( _sLeftImage.c_str(), &imageL );
    image::ReadImage( _sRightImage.c_str(), &imageR );

    image::Image<float> imgA ( imageL.GetMat().cast<float>() );
    image::Image<float> imgB(imageR.GetMat().cast<float>());

    std::vector<Pair> matchesFiltered, matchesPair;

    for (const auto & iter_match : _vec_PutativeMatches)
    {
      matchesPair.push_back( {iter_match.i_, iter_match.j_} );
    }

    std::vector<double> vec_score;

    //In order to illustrate the gvld(or vld)-consistant neighbors, the following two parameters has been externalized as inputs of the function KVLD.
    openMVG::Mat E = openMVG::Mat::Ones( _vec_PutativeMatches.size(), _vec_PutativeMatches.size() ) * ( -1 );
    // gvld-consistancy matrix, intitialized to -1,  >0 consistancy value, -1=unknow, -2=false
    std::vector<bool> valide( _vec_PutativeMatches.size(), true );// indices of match in the initial matches, if true at the end of KVLD, a match is kept.

    size_t it_num = 0;
    KvldParameters kvldparameters;//initial parameters of KVLD
    //kvldparameters.K = 5;
    while (
      it_num < 5 &&
      kvldparameters.inlierRate >
      KVLD(
        imgA, imgB,
        _vec_featsL, _vec_featsR,
        matchesPair, matchesFiltered,
        vec_score, E, valide, kvldparameters ) )
    {
      kvldparameters.inlierRate /= 2;
      std::cout<<"low inlier rate, re-select matches with new rate="<<kvldparameters.inlierRate<<std::endl;
      kvldparameters.K = 2;
      it_num++;
    }

    bool bOk = false;
    if (!matchesPair.empty())
    {
      // Get mask
      getKVLDMask(
        &maskLeft, &maskRight,
        _vec_featsL, _vec_featsR,
        matchesPair,
        valide,
        E);
      bOk = true;
    }
    else
    {
      maskLeft.fill( 0 );
      maskRight.fill( 0 );
    }

    return bOk;
  }

private:
  // Left and Right features
  std::vector<features::SIOPointFeature > _vec_featsL, _vec_featsR;
  // Left and Right corresponding index (putatives matches)
  std::vector<matching::IndMatch> _vec_PutativeMatches;
};

}  // namespace color_harmonization
}  // namespace openMVG

#endif  // OPENMVG_COLOR_HARMONIZATION_SELECTION_VLDSEGMENT_HPP
