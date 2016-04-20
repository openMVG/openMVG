
// Copyright (c) 2014 openMVG authors.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_COLORHARMONIZATION_MATCHESPOINTS_H
#define OPENMVG_COLORHARMONIZATION_MATCHESPOINTS_H

#include "openMVG/color_harmonization/selection_interface.hpp"
#include "openMVG/matching/indMatch.hpp"
#include "openMVG/features/features.hpp"

#include <vector>

namespace openMVG {
namespace color_harmonization {

class commonDataByPair_MatchedPoints  : public commonDataByPair
{
public:
  commonDataByPair_MatchedPoints(const std::string & sLeftImage,
                                 const std::string & sRightImage,
                                 const std::vector< matching::IndMatch >& vec_PutativeMatches,
                                 const std::vector< features::SIOPointFeature >& vec_featsL,
                                 const std::vector< features::SIOPointFeature >& vec_featsR,
                                 const size_t radius = 1 ):
     commonDataByPair( sLeftImage, sRightImage ),
     _vec_PutativeMatches( vec_PutativeMatches ),
     _vec_featsL( vec_featsL ), _vec_featsR( vec_featsR ), _radius( radius )
  {}

  ~commonDataByPair_MatchedPoints() override = default ; 

  /**
   * Fill mask from corresponding points (each point pictured by a disk of radius _radius)
   *
   * \param[out] maskLeft Mask of the left image (initialized to corresponding image size).
   * \param[out] maskRight  Mask of the right image  (initialized to corresponding image size).
   *
   * \return True if some pixel have been set to true.
   */
  bool computeMask( image::Image< unsigned char > & maskLeft, image::Image< unsigned char > & maskRight ) override 
  {
    maskLeft.fill(0);
    maskRight.fill(0);
    for( std::vector< matching::IndMatch >::const_iterator
          iter_putativeMatches = _vec_PutativeMatches.begin();
          iter_putativeMatches != _vec_PutativeMatches.end();
          ++iter_putativeMatches )
    {
      const features::SIOPointFeature & L = _vec_featsL[ iter_putativeMatches->i_ ];
      const features::SIOPointFeature & R = _vec_featsR[ iter_putativeMatches->j_ ];

      image::FilledCircle( L.x(), L.y(), ( int )_radius, ( unsigned char ) 255, &maskLeft );
      image::FilledCircle( R.x(), R.y(), ( int )_radius, ( unsigned char ) 255, &maskRight );
    }
    return _vec_PutativeMatches.size() > 0;
  }

private:
  size_t _radius;
  std::vector< matching::IndMatch > _vec_PutativeMatches;
  std::vector< features::SIOPointFeature > _vec_featsL;
  std::vector< features::SIOPointFeature > _vec_featsR;
};

}  // namespace color_harmonization
}  // namespace openMVG

#endif  // OPENMVG_COLORHARMONIZATION_MATCHESPOINTS_H
