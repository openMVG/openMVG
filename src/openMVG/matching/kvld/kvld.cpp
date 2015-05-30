/** @Main KVLD algorithm implementation
 ** @Containing scale image pyramid, VLD structure and KVLD algorithm
 ** @author Zhe Liu
 **/

/*
Copyright (C) 2011-12 Zhe Liu and Pierre Moulon.
All rights reserved.

This file is part of the KVLD library and is made available under
the terms of the BSD license ( see the COPYING file).
*/

#include "kvld.h"
#include "algorithm.h"
#include <functional>
#include <numeric>
#include <openMVG/image/image.hpp>

using namespace std;
using namespace openMVG;
using namespace openMVG::image;

ImageScale::ImageScale( const Image< float >& I, double r )
{
  IntegralImages inter( I );
  radius_size = r;
  step = sqrt( 2.0 );
  int size = max( I.Width(),I.Height() );

  int number= int( log( size / r ) / log( 2.0 ) ) + 1;
  angles.resize( number );
  magnitudes.resize( number );
  ratios.resize( number );

  GradAndNorm( I, angles[ 0 ], magnitudes[ 0 ] );
  ratios[ 0 ] = 1;

#ifdef OPENMVG_USE_OPENMP
#pragma omp parallel for
#endif
  for( int k = 1; k < number; k++ )
  {
    Image< float > I2;
    double ratio = 1 * pow( step, k );
    I2.resize( int( I.Width() / ratio ), int( I.Height() / ratio ) );
    angles[ k ].resize( int( I.Width() / ratio ), int( I.Height() / ratio ) );
    magnitudes[ k ].resize( int( I.Width() / ratio ), int( I.Height() / ratio ) );

    for( int i = 0; i < I2.Width(); i++ )
    {
      for( int j = 0; j < I2.Height(); j++ )
      {
        I2( j, i ) = inter( double( i + 0.5 ) * ratio, double( j + 0.5 ) * ratio, ratio );
      }
    }
    GradAndNorm( I2,angles[ k ], magnitudes[ k ] );
    ratios[ k ] = ratio;
  }
}

void ImageScale::GradAndNorm( const Image< float >& I, Image< float >& angle, Image< float >& m )
{
  angle = Image< float >( I.Width(), I.Height() );
  m = Image< float >( I.Width(), I.Height() );
  angle.fill( 0 );
  m.fill( 0 );
#ifdef OPENMVG_USE_OPENMP
#pragma omp parallel for
#endif
  for( int y = 1; y < I.Height() - 1; y++ )
  {
  for( int x = 1; x < I.Width() - 1; x++ )
  {
      const float gx = I( y, x + 1 ) - I( y, x - 1 );
      const float gy = I( y + 1, x ) - I( y - 1, x );

      if( !anglefrom( gx, gy, angle( y, x ) ) )
        angle( y, x ) = -1;
      m( y, x ) = sqrt( gx * gx + gy * gy );
    }
  }
}

int ImageScale::getIndex( const double r )const
{
  const double step = sqrt( 2.0 );

  if( r <= radius_size ) return 0;
  else
  {
    double range_low = radius_size;
    int index = 0;
    while( r > range_low * step )
    {
      ++index;
      range_low *= step;
    }
    return std::min(int(angles.size()-1), index);
  }
}

template< typename T >
VLD::VLD( const ImageScale& series, T const& P1, T const& P2 ) : contrast( 0.0 )
{
  //============== initializing============//
  principleAngle.fill( 0 );
  descriptor.fill( 0 );
  weight.fill( 0 );

  begin_point[ 0 ] = P1.x();
  begin_point[ 1 ] = P1.y();
  end_point[ 0 ]   = P2.x();
  end_point[ 1 ]   = P2.y();

  const float dy = float( end_point[ 1 ] - begin_point[ 1 ] );
  const float dx = float( end_point[ 0 ] - begin_point[ 0 ] );
  distance = sqrt( dy * dy + dx * dx );

  if( distance == 0 )
    cerr<<"Two SIFT points have the same coordinate"<<endl;

  const float radius = max( distance / float( dimension + 1 ), 2.0f );//at least 2

  const double mainAngle = get_orientation();//absolute angle

  const int image_index = series.getIndex( radius );

  const Image< float > & ang = series.angles[ image_index ];
  const Image< float > & m   = series.magnitudes[ image_index ];
  const double ratio = series.ratios[ image_index ];

  const int w = m.Width();
  const int h = m.Height();
  const float r = float( radius / ratio );
  const float sigma2 = r * r;
  //======calculating the descriptor=====//

  double statistic[ binNum ];
  for( int i = 0; i < dimension; i++ )
  {
    fill_n( statistic, binNum, 0.0);

    float xi = float( begin_point[ 0 ] + float( i + 1 ) / ( dimension + 1 ) * ( dx ) );
    float yi = float( begin_point[ 1 ] + float( i + 1 ) / ( dimension + 1 ) * ( dy ) );
    xi /= float( ratio );
    yi /= float( ratio );

    for( int y = int( yi - r ); y <= int( yi + r + 0.5 ); y++ )
    {
      for( int x = int( xi - r ); x <= int( xi + r + 0.5 ); x++ )
      {
        float d = point_distance( xi, yi, float( x ), float( y ) );
        if( d <= r && inside( w, h, x, y, 1 ) )
        {
          //================angle and magnitude==========================//
          double angle;
          if( ang( y, x ) >= 0 )
            angle = ang( y, x ) - mainAngle;//relative angle
          else angle = 0.0;

          //cout<<angle<<endl;
          while( angle < 0 )
            angle += 2 * PI_;
          while( angle >= 2 * PI_)
            angle -= 2 * PI_;

          //===============principle angle==============================//
          const int index = int( angle * binNum / ( 2 * PI_ ) + 0.5 );

          double Gweight = exp( -d * d / 4.5 / sigma2 ) * ( m( y, x ) );
          if( index < binNum )
            statistic[ index ] += Gweight;
          else // possible since the 0.5
            statistic[ 0 ] += Gweight;

          //==============the descriptor===============================//
          const int index2 = int( angle * subdirection / ( 2 * PI_ ) + 0.5 );
          assert( index2 >= 0 && index2 <= subdirection );

          if( index2 < subdirection )
            descriptor[ subdirection * i + index2 ] += Gweight;
          else descriptor[ subdirection * i ] += Gweight;// possible since the 0.5
        }
      }
    }
    //=====================find the biggest angle of ith SIFT==================//
    int index;
    int second_index;
    max( statistic, weight[ i ], binNum, index, second_index );
    principleAngle[ i ] = index;
  }

  normalize_weight( descriptor );

  contrast  = weight.array().sum();
  contrast /= distance / ratio;
  normalize_weight( weight );
}

float KVLD( const Image< float >& I1,
            const Image< float >& I2,
            const std::vector<features::SIOPointFeature> & F1,
            const std::vector<features::SIOPointFeature> & F2,
            const vector< Pair >& matches,
            vector< Pair >& matchesFiltered,
            vector< double >& score,
            openMVG::Mat& E,
            vector< bool >& valide,
            KvldParameters& kvldParameters )
{
  matchesFiltered.clear();
  score.clear();
  ImageScale Chaine1( I1 );
  ImageScale Chaine2( I2 );

  cout << "Image scale-space complete..." << endl;

  const float range1 = getRange( I1, min( F1.size(), matches.size() ), kvldParameters.inlierRate );
  const float range2 = getRange( I2, min( F2.size(), matches.size() ), kvldParameters.inlierRate );

  const size_t size = matches.size();

  //================distance map construction, foruse of selecting neighbors===============//
  cout << "computing distance maps" << endl;

  bool bPrecomputedDist = false;

  openMVG::Matf dist1, dist2;
  if( bPrecomputedDist )
  {
    dist1 = openMVG::Matf::Zero( F1.size(), F1.size() );
    dist2 = openMVG::Matf::Zero( F2.size(), F2.size() );

    for( int a1 = 0; a1 < F1.size(); ++a1 )
      for( int a2 = a1; a2 < F1.size(); ++a2 )
        dist1( a1, a2 ) = dist1( a2, a1 ) = point_distance( F1[ a1 ], F1[ a2 ] );

    for( int b1 = 0; b1 < F2.size(); ++b1 )
      for( int b2 = b1; b2 < F2.size(); ++b2 )
        dist2( b1, b2 ) = dist2( b2, b1 ) = point_distance( F2[ b1 ], F2[ b2 ] );
  }

  fill( valide.begin(), valide.end(), true );
  vector< double > scoretable( size, 0.0 );
  vector< size_t > result( size, 0 );


//============main iteration formatch verification==========//
//    cout<<"main iteration";
  bool change = true;

  while( change )
  {
    change = false;

    fill( scoretable.begin(), scoretable.end(), 0.0 );
    fill( result.begin(), result.end(), 0 );
    //========substep 1: search foreach match its neighbors and verify if they are gvld-consistent ============//
    for( int it1 = 0; it1 < size - 1; it1++ )
    {
      if( valide[ it1 ] )
      {
        size_t a1 = matches[ it1 ].first, b1 = matches[ it1 ].second;

        for( int it2 = it1 + 1; it2 < size; it2++ )
          if(valide[ it2 ])
          {
            size_t a2 = matches[ it2 ].first, b2 = matches[ it2 ].second;

            bool bOk = false;
            if( bPrecomputedDist )
              bOk = ( dist1( a1, a2 ) > min_dist && dist2( b1, b2 ) > min_dist
                 && ( dist1( a1, a2 ) < range1   || dist2( b1, b2 ) < range2 ) );
            else
              bOk = ( point_distance( F1[ a1 ], F1[ a2 ] ) > min_dist && point_distance( F2[ b1 ], F2[ b2 ] ) > min_dist &&
                    ( point_distance( F1[ a1 ], F1[ a2 ] ) < range1   || point_distance( F2[ b1 ], F2[ b2 ] ) < range2 ) );
            if( bOk )
            {
              if( E( it1, it2 ) == -1 )
              { //update E ifunknow
                E( it1, it2 ) = -2;
                E( it2, it1 ) = -2;

                if( !kvldParameters.geometry || consistent( F1[ a1 ], F1[ a2 ], F2[ b1 ], F2[ b2 ] ) < distance_thres )
                {
                  VLD vld1( Chaine1, F1[ a1 ], F1[ a2 ] );
                  VLD vld2( Chaine2, F2[ b1 ], F2[ b2 ] );
                  //vld1.test();
                  double error = vld1.difference( vld2 );
                  //cout<<endl<<it1<<" "<<it2<<" "<<dist1(a1,a2)<<" "<< dist2(b1,b2)<<" "<<error<<endl;
                  if( error < juge )
                  {
                    E( it1, it2 ) = ( float ) error;
                    E( it2, it1 ) = ( float ) error;
                    //cout<<E(it2,it1)<<endl;
                  }
                }
              }
              if( E( it1, it2 ) >= 0 )
              {
                result[ it1 ] += 1;
                result[ it2 ] += 1;
                scoretable[ it1 ] += double( E( it1, it2 ) );
                scoretable[ it2 ] += double( E( it1, it2 ) );
                if( result[ it1 ] >= max_connection )
                  break;
              }
            }
          }
      }
    }

    //========substep 2: remove false matches by K gvld-consistency criteria ============//
    for( int it = 0; it < size; it++ )
    {
      if( valide[ it ] && result[ it ] < kvldParameters.K )
      {
        valide[ it ] = false;
        change = true;
      }
    }
    //========substep 3: remove multiple matches to a same point by keeping the one with the best average gvld-consistency score ============//
    if( uniqueMatch )
      for( int it1 = 0; it1 < size - 1; it1++ )
        if( valide[ it1 ]) {
          size_t a1 = matches[ it1 ].first;
          size_t b1 = matches[ it1 ].second;

          for( int it2 = it1 + 1; it2 < size; it2++ )
            if( valide[ it2 ] )
            {
              size_t a2 = matches[ it2 ].first;
              size_t b2 = matches[ it2 ].second;

              if( a1 == a2 || b1 == b2
                  || ( F1[ a1 ].x() == F1[ a2 ].x() && F1[ a1 ].y() == F1[ a2 ].y() &&
                     ( F2[ b1 ].x() != F2[ b2 ].x() || F2[ b1 ].y() != F2[ b2 ].y() ) )
                  || ( ( F1[ a1 ].x() != F1[ a2 ].x() || F1[ a1 ].y() != F1[ a2 ].y() ) &&
                         F2[ b1 ].x() == F2[ b2 ].x() && F2[ b1 ].y() == F2[ b2 ].y() ) )
              {
                //cardinal comparison
                if( result[ it1 ] > result[ it2 ] )
                {
                  valide[ it2 ] = false;
                  change = true;
                }
                else if( result[ it1 ] < result[ it2 ] )
                {
                  valide[ it1 ] = false;
                  change = true;
                }
                else if( result[ it1 ] == result[ it2 ] )
                {
                  //score comparison
                  if( scoretable[ it1 ] > scoretable[ it2 ] )
                  {
                    valide[ it1 ] = false;
                    change = true;
                  }
                  else if( scoretable[ it1 ] < scoretable[ it2 ] )
                  {
                    valide[ it2 ] = false;
                    change = true;
                  }
                }
              }
            }
        }
    //========substep 4: ifgeometric verification is set, re-score matches by geometric-consistency, and remove poorly scored ones ============================//
    if( uniqueMatch && kvldParameters.geometry )
    {
      for( int i = 0; i < size; i++ )
        scoretable[ i ]=0;

      vector< bool > switching;
      for( int i = 0; i < size; i++ )
        switching.push_back( false );

      for( int it1 = 0; it1 < size; it1++ )
      {
        if( valide[ it1 ] )
        {
          size_t a1 = matches[ it1 ].first, b1 = matches[ it1 ].second;
          float index = 0.0f;
          int good_index = 0;
          for( int it2 = 0; it2 < size; it2++ )
          {
            if( it1 != it2 && valide[ it2 ] )
            {
              size_t a2 = matches[ it2 ].first;
              size_t b2 = matches[ it2 ].second;

              bool bOk = false;
              if( bPrecomputedDist )
                bOk = ( dist1( a1, a2 ) > min_dist && dist2( b1, b2 ) > min_dist &&
                      ( dist1( a1, a2 ) < range1   || dist2( b1, b2 ) < range2 ) );
              else
                bOk = ( point_distance( F1[ a1 ], F1[ a2 ] ) > min_dist && point_distance( F2[ b1 ], F2[ b2 ] ) > min_dist &&
                      ( point_distance( F1[ a1 ], F1[ a2 ] ) < range1   || point_distance( F2[ b1 ], F2[ b2 ] ) < range2 ) );
              if( bOk )
              {
                float d = consistent( F1[ a1 ], F1[ a2 ], F2[ b1 ], F2[ b2 ] );
                scoretable[ it1 ] += d;
                index += 1;
                if( d < distance_thres )
                  good_index++;
              }
            }
          }
          scoretable[ it1 ] /= index;
          if( good_index < 0.3f * float( index ) && scoretable[ it1 ] > 1.2 )
          {
            switching[ it1 ] = true;
            change = true;
          }
        }
      }
      for( int it1 = 0; it1 < size; it1++ )
        if( switching[ it1 ] )
          valide[ it1 ] = false;
    }
  }
  //=============== generating output list ===================//
  for( int it = 0; it < size; it++ )
    if( valide[ it ] )
    {
      matchesFiltered.push_back( matches[ it ] );
      score.push_back( scoretable[ it ] );
    }
  return float( matchesFiltered.size() ) / matches.size();
}

