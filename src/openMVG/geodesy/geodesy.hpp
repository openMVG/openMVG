// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2016 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_GEODESY_GEODESY_HPP
#define OPENMVG_GEODESY_GEODESY_HPP

#include "openMVG/numeric/numeric.h"

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
 ** Convert WGS84 lon,lat,alt data to UTM data (Universal Transverse Mercator)
 ** @param lat Latitude in degree
 ** @param lon Longitude in degree
 ** @param alt Altitude relative to the WGS84 ellipsoid
 ** @return UTM corresponding coordinates
 **/
Vec3 lla_to_utm
(
   double lat,
   double lon,
   double alt,
   double a = WGS84_A,
   double b = WGS84_B
)
{
  a /= 1000; // meters to kilometers
  b /= 1000; // meters to kilometers

  /// CONSTANTS
  static const double N0_n = 0;
  static const double N0_s = 1e4;
  static const double E0 = 5e2;
  static const double k0 = 0.9996;
  const double f = (a - b) / a;

  const double n    = f / (2 - f);
  const double n_2  = n   * n;
  const double n_3  = n_2 * n;
  const double n_4  = n_3 * n;
  const double n_5  = n_4 * n;
  const double n_6  = n_5 * n;
  const double n_7  = n_6 * n;
  const double n_8  = n_7 * n;
  const double n_9  = n_8 * n;
  const double n_10 = n_9 * n;

  const int lon_zone = 1 + floor((lon + 180) / 6);

  const double lon_0 = D2R(3 + 6 * (lon_zone - 1) - 180);

  lat = D2R(lat);
  lon = D2R(lon);

  const double A = a / (1 + n) * (1 + n_2/4 + n_4/64 + n_6/256 + n_8*25.0/16384.0 + n_10*49.0/65536.0);

  const double a1 = (1.0/2.0)*n - (2.0/3.0)*n_2 + (5.0/16.0)*n_3 + (41.0/180.0)*n_4 - (127.0/288.0)*n_5 + (7891.0/37800.0)*n_6 + (72161.0/387072.0)*n_7 - (18975107.0/50803200.0)*n_8 + (60193001.0/290304000.0)*n_9 + (134592031.0/1026432000.0)*n_10;
  const double a2 = (13.0/48.0)*n_2 - (3.0/5.0)*n_3 + (557.0/1440.0)*n_4 + (281.0/630.0)*n_5 - (1983433.0/1935360.0)*n_6 + (13769.0/28800.0)*n_7 + (148003883.0/174182400.0)*n_8 - (705286231.0/465696000.0)*n_9 + (1703267974087.0/3218890752000.0)*n_10;
  const double a3 = (61.0/240.0)*n_3 - (103.0/140.0)*n_4 + (15061.0/26880.0)*n_5 + (167603.0/181440.0)*n_6 - (67102379.0/29030400.0)*n_7 + (79682431.0/79833600.0)*n_8 + (6304945039.0/2128896000.0)*n_9 - (6601904925257.0/1307674368000.0)*n_10;
  const double a4 = (49561.0/161280.0)*n_4 - (179.0/168.0)*n_5 + (6601661.0/7257600.0)*n_6 + (97445.0/49896.0)*n_7 - (40176129013.0/7664025600.0)*n_8 + (138471097.0/66528000.0)*n_9 + (48087451385201.0/5230697472000.0)*n_10;
  const double a5 = (34729.0/80640.0)*n_5 - (3418889.0/1995840.0)*n_6 + (14644087.0/9123840.0)*n_7 + (2605413599.0/622702080.0)*n_8 - (31015475399.0/2583060480.0)*n_9 + (5820486440369.0/1307674368000.0)*n_10;
  const double a6 = (212378941.0/319334400.0)*n_6 - (30705481.0/10378368.0)*n_7 + (175214326799.0/58118860800.0)*n_8 + (870492877.0/96096000.0)*n_9 - (1328004581729009.0/47823519744000.0)*n_10;
  const double a7 = (1522256789.0/1383782400.0)*n_7 - (16759934899.0/3113510400.0)*n_8 + (1315149374443.0/221405184000.0)*n_9 + (71809987837451.0/3629463552000.0)*n_10;
  const double a8 = (1424729850961.0/743921418240.0)*n_8 - (256783708069.0/25204608000.0)*n_9 + (2468749292989891.0/203249958912000.0)*n_10;
  const double a9 = (21091646195357.0/6080126976000.0)*n_9 - (67196182138355857.0/3379030566912000.0)*n_10;
  const double a10 = (77911515623232821.0/12014330904576000.0)*n_10;

  const double t = sinh(atanh(sin(lat)) - 2*sqrt(n)/(1+n) * atanh(2*sqrt(n)/(1+n)*sin(lat)));
  const double xi = atan(t/cos(lon-lon_0));
  const double eta = atanh(sin(lon-lon_0) / sqrt(1+t*t));

  const double N0 = (lat > 0 ? N0_n : N0_s);

  const double E = E0 + k0 * A * (eta + a1*cos(2*1*xi)*sinh(2*1*eta) + a2*cos(2*2*xi)*sinh(2*2*eta) + a3*cos(2*3*xi)*sinh(2*3*eta) + a4*cos(2*4*xi)*sinh(2*4*eta) + a5*cos(2*5*xi)*sinh(2*5*eta) + a6*cos(2*6*xi)*sinh(2*6*eta) + a7*cos(2*7*xi)*sinh(2*7*eta) + a8*cos(2*8*xi)*sinh(2*8*eta) + a9*cos(2*9*xi)*sinh(2*9*eta) + a10*cos(2*10*xi)*sinh(2*10*eta));
  const double N = N0 + k0 * A * (xi + a1*sin(2*1*xi)*cosh(2*1*eta) + a2*sin(2*2*xi)*cosh(2*2*eta) + a3*sin(2*3*xi)*cosh(2*3*eta) + a4*sin(2*4*xi)*cosh(2*4*eta) + a5*sin(2*5*xi)*cosh(2*5*eta) + a6*sin(2*6*xi)*cosh(2*6*eta) + a7*sin(2*7*xi)*cosh(2*7*eta) + a8*sin(2*8*xi)*cosh(2*8*eta) + a9*sin(2*9*xi)*cosh(2*9*eta) + a10*sin(2*10*xi)*cosh(2*10*eta));

  // Scale E,N from kilometers to meters
  return Vec3(E * 1000, N * 1000, alt);
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


#endif // OPENMVG_GEODESY_GEODESY_HPP
