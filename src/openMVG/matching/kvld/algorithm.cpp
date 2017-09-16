/** @basic structures implementation
 ** @author Zhe Liu
 **/

/*
Copyright (C) 2011-12 Zhe Liu and Pierre Moulon.
All rights reserved.

This file is part of the KVLD library and is made available under
the terms of the BSD license (see the COPYING file).
*/

#include "algorithm.h"

#include "openMVG/features/feature.hpp"

IntegralImages::IntegralImages(const openMVG::image::Image< float >& I)
{
  map.resize( I.Width() + 1, I.Height() + 1 );
  map.fill( 0 );
  for (int y = 0; y < I.Height(); y++ )
  {
    for (int x = 0; x < I.Width(); x++ )
    {
      map( y + 1, x + 1 ) = double( I( y, x ) ) + map( y, x + 1 ) + map( y + 1, x ) - map( y, x );
    }
  }
}

float getRange(
  const openMVG::image::Image< float >& I,
  int a,
  const float p )
{
  float range = sqrt( float( 3.f * I.Height() * I.Width() ) / ( p * a * PI_ ) );
  return range;
}


//=============================IO interface======================//

std::ofstream& writeDetector( std::ofstream& out, const openMVG::features::SIOPointFeature& feature )
{
  out << feature.x() << " "
    << feature.y() << " "
    << feature.scale() << " "
    << feature.orientation() <<std::endl;
  return out;
}

std::ifstream& readDetector( std::ifstream& in, openMVG::features::SIOPointFeature& point )
{
  in >> point.x()
    >> point.y()
    >> point.scale()
    >> point.orientation();
  return in;
}
