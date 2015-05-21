#ifndef KVLD_ALGORITHM_H
#define KVLD_ALGORITHM_H

/** @basic structures implementation
 ** @author Zhe Liu
 **/

/*
Copyright (C) 2011-12 Zhe Liu and Pierre Moulon.
All rights reserved.

This file is part of the KVLD library and is made available under
the terms of the BSD license (see the COPYING file).
*/

#include <openMVG/numeric/numeric.h>
#include <openMVG/image/image_container.hpp>
#include <openMVG/matching/indMatch.hpp>
#include <openMVG/features/feature.hpp>
#include <openMVG/types.hpp>

#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>
#include <numeric>
#include <memory>
#include <algorithm>
#include <functional>


const float PI_ = 4.0 * atan( 1.0f );

//============================== simplified structure of a point=============================//
//if you set KVLD geometry verification to false, you only need to fill x and y in a point structure
struct PointS
{
	float x, y; // point position
	float scale; // point scale
	float angle; // point orientation

	PointS( float x = (0.f), float y = (0.f)):
      x( x ), y( y ), scale(0.f), angle(0.f){}
	PointS( const float& x, const float& y,const float& angle,const float& scale):
      x( x ), y( y ), angle( angle ), scale( scale ){}
};

//===================================== integral image ====================================//
//It is used to efficiently construct the pyramid of scale images in KVLD
struct IntegralImages
{
  openMVG::image::Image< double > map;

  IntegralImages(const openMVG::image::Image< float >& I);

  inline double operator()( double x1, double y1, double x2, double y2 )const
  {
		return get( x2, y2 ) - get( x1, y2 ) - get( x2, y1 ) + get( x1, y1 );
	}
  inline double operator()( double x, double y, double size ) const
  {
    double window = 0.5 * size;
    return ( get( x + window, y + window ) - get( x - window, y + window ) - get( x + window, y - window ) + get( x - window, y - window ) ) / ( 4 * window * window );
  }
private :
  inline double get( double x, double y )const
  {
		int ix = int( x );
		int iy = int( y );
		double dx = x - ix;
		double dy = y - iy;
		if( dx == 0 && dy == 0 )
			return map( iy, ix );
		if( dx == 0 )
			return map( iy, ix ) * ( 1 - dy ) + map( iy + 1, ix ) * dy;
		if( dy == 0 )
			return map( iy, ix ) * ( 1 - dx ) + map( iy, ix + 1 ) * dx;

		return map( iy, ix ) * ( 1 - dx ) * ( 1 - dy ) +
            map( iy + 1, ix ) * dy * ( 1 - dx ) +
            map( iy, ix + 1 ) * ( 1 - dy ) * dx +
            map( iy + 1, ix + 1 ) * dx * dy;
	}
};

//=============================IO interface ======================//

std::ofstream& writeDetector( std::ofstream& out, const openMVG::features::SIOPointFeature& vect );
std::ifstream& readDetector( std::ifstream& in, openMVG::features::SIOPointFeature& point );
//======================================elemetuary operations================================//
template < typename T >
inline T point_distance( const T x1, const T y1, const T x2, const T y2 )
{//distance of points
	float a = x1 - x2;
	float b = y1 - y2;
	return sqrt( a * a + b * b );
}

template < typename T >
inline float point_distance( const T& P1, const T& P2 )
{//distance of points
		 return point_distance< float >( P1.x(), P1.y(), P2.x(), P2.y() );
}

inline bool inside( int w, int h, int x,int y, double radius )
{
	return( x - radius >= 0 && y - radius >= 0 && x + radius < w && y + radius < h );
}

inline bool anglefrom( const float& x, const float& y, float& angle )
{
	if( x != 0 )
		angle = atan( y / x );
	else if( y > 0 )
		angle = PI_ / 2;
	else if( y < 0 )
		angle =- PI_ / 2;
	else return false;

	if( x < 0 )
		angle += PI_;
	while( angle < 0 )
		angle += 2 * PI_;
  while( angle >= 2 * PI_ )
		angle -= 2 * PI_;
	assert( angle >= 0 && angle < 2 * PI_ );
	return true;
}

inline double angle_difference( const double angle1, const double angle2 )
{
	double angle = angle1 - angle2;
	while( angle <  0 ) angle += 2 * PI_;
	while( angle >= 2 * PI_ )	angle -= 2 * PI_;

	assert(angle <= 2 * PI_ && angle >= 0 );
	return std::min( angle, 2 * PI_ - angle );
}

inline void max( double* list,double& weight, int size, int& index, int& second_index )
{
	index = 0;
	second_index = -1;
	double best = list[ index ] - list[ index + size / 2 ];

	for( int i = 0; i < size; i++ )
  {
			double value;
			if( i < size / 2) value = list[ i ] - list[ i + size / 2 ];
			else value = list[ i ] - list[ i - size / 2 ];

			if( value > best )
      {
				best = value;
				second_index = index;
				index = i;
			}
	}
	weight = best;
}

template< typename ARRAY >
inline void normalize_weight( ARRAY & weight )
{
  double total = weight.array().sum();
	if( total != 0 )
  	for( int i = 0; i < weight.size(); i++ )
	  	weight[ i ] /= total;
}

template< typename T >
inline float consistent( const T& a1, const T& a2, const T& b1, const T& b2 )
{
	float ax = float( a1.x() - a2.x() );
	float ay = float( a1.y() - a2.y() );
	float bx = float( b1.x() - b2.x() );
	float by = float( b1.y() - b2.y() );

	float angle1 = float( b1.orientation() - a1.orientation() );
	float angle2 = float( b2.orientation() - a2.orientation() );

	float ax1 = cos( angle1 ) * ax - sin( angle1 ) * ay;
	ax1 *= float( b1.scale() / a1.scale() );
	float ay1 = sin( angle1 ) * ax + cos( angle1 ) * ay;
	ay1 *= float( b1.scale() / a1.scale() );
	float d1 = sqrt( ax1 * ax1 + ay1 * ay1 );
	float d1_error = sqrt( ( ax1 - bx ) * ( ax1 - bx ) + ( ay1 - by ) * ( ay1 - by ) );

	float ax2 = float( cos( angle2 ) * ax - sin( angle2 ) * ay );
	ax2 *= float( b2.scale() / a2.scale() );

	float ay2 = float( sin( angle2 ) * ax + cos( angle2 ) * ay );
	ay2 *= float( b2.scale() / a2.scale() );
	float d2 = sqrt( ax2 * ax2 + ay2 * ay2 );
	float d2_error = sqrt( ( ax2 - bx ) * ( ax2 - bx ) + ( ay2 - by ) * ( ay2 - by ) );
	float d = std::min( d1_error / std::min( d1, point_distance( b1, b2 ) ), d2_error / std::min( d2, point_distance( b1, b2 ) ) );
	return d;
}
float getRange(const openMVG::image::Image< float >& I, int a, const float p);

#endif //KVLD_ALGORITHM_H
