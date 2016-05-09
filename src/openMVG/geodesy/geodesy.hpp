// Copyright (c) 2016 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <openMVG/numeric/numeric.h>

namespace openMVG
{
/**
* @brief namespace containing various classes and functions used to deal with geodesy transformation
*/
namespace geodesy
{

// WGS84 Ellipsoid
static const double WGS84_A = 6378137.0;      // major axis
static const double WGS84_B = 6356752.314245; // minor axis
static const double WGS84_E = 0.0818191908;   // first eccentricity

/**
 ** Convert WGS84 lon,lat,alt data to ECEF data (Earth Centered Earth Fixed)
 ** @param lat Latitude in degree
 ** @param lon Longitude in degree
 ** @param alt Altitude relative to the WGS84 ellipsoid
 ** @return ECEF corresponding coordinates
 **/
Vec3 lla_to_ecef
(
  double lat,
  double lon,
  double alt
)
{
  const double clat = cos( D2R(lat) );
  const double slat = sin( D2R(lat) );
  const double clon = cos( D2R(lon) );
  const double slon = sin( D2R(lon) );

  const double a2 = Square(WGS84_A);
  const double b2 = Square(WGS84_B);

  const double L = 1.0 / sqrt(a2 * Square(clat) + b2 * Square(slat));
  const double x = (a2 * L + alt) * clat * clon;
  const double y = (a2 * L + alt) * clat * slon;
  const double z = (b2 * L + alt) * slat;

  return Vec3(x, y, z);
}

/**
 ** Convert ECEF (XYZ) to lon,lat,alt values for the WGS84 ellipsoid
 ** @param x X ECEF coordinate
 ** @param y Y ECEF coordinate
 ** @param z Z ECEF coordinate
 ** @return LLA corresponding coordinates
 **/ 
// http://fr.mathworks.com/matlabcentral/newsreader/view_thread/142629
Vec3 ecef_to_lla
(
  double x,
  double y,
  double z
)
{
  const double b = sqrt(WGS84_A*WGS84_A*(1-WGS84_E*WGS84_E));
  const double ep = sqrt((WGS84_A*WGS84_A-b*b)/(b*b));
  const double p = sqrt(x*x+y*y);
  const double th = atan2(WGS84_A*z,b*p);
  const double lon = atan2(y,x);
  const double lat = atan2((z+ep*ep*b* pow(sin(th),3)),(p-WGS84_E*WGS84_E*WGS84_A*pow(cos(th),3)));
  const double N = WGS84_A/sqrt(1-WGS84_E*WGS84_E*sin(lat)*sin(lat));
  const double alt = p/cos(lat)-N;

  return Vec3(R2D(lat), R2D(lon), alt);
}

} // namespace geodesy
} // namespace openMVG

