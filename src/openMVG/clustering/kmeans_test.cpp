// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Romuald Perrot.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/clustering/kmeans.hpp"

#include "testing/testing.h"

#include <random>
#include <vector>

using namespace openMVG ;
using namespace clustering ;

TEST( clustering , threeClustersVec2 )
{
  std::vector< Vec2 > pts ;

  const Vec2 center1( -5.0 , -5.0 ) ;
  const Vec2 center2( 5.0 , 5.0 ) ;
  const Vec2 center3( 0.0 , 0.0 ) ;

  const int nb_pts_1 = 10000 ;
  const int nb_pts_2 = 10000 ;
  const int nb_pts_3 = 10000 ;

  std::mt19937_64 rng ;
  std::uniform_real_distribution<double> distrib( -2.0 , 2.0 ) ;

  // Generate points on each centers
  for( int id_pt_1 = 0 ; id_pt_1 < nb_pts_1 ; ++id_pt_1 )
  {
    const double dx = distrib( rng ) ;
    const double dy = distrib( rng ) ;

    pts.push_back( center1 + Vec2( dx , dy ) ) ;
  }
  for( int id_pt_2 = 0 ; id_pt_2 < nb_pts_2 ; ++id_pt_2 )
  {
    const double dx = distrib( rng ) ;
    const double dy = distrib( rng ) ;

    pts.push_back( center2 + Vec2( dx , dy  ) ) ;
  }
  for( int id_pt_3 = 0 ; id_pt_3 < nb_pts_3 ; ++id_pt_3 )
  {
    const double dx = distrib( rng ) ;
    const double dy = distrib( rng ) ;

    pts.push_back( center3 + Vec2( dx , dy  ) ) ;
  }

  std::vector< uint32_t > ids ;
  std::vector< Vec2 > centers ;
  KMeans( pts , ids , centers , 3 ) ;

  std::cout << "Centers : " << std::endl ;
  for( const auto & it : centers )
  {
    std::cout << it << std::endl ;
  }

  // First set
  const uint32_t id_0 = ids[ 0 ] ;
  for( int id_pt_1 = 1 ; id_pt_1 < nb_pts_1 ; ++id_pt_1 )
  {
    EXPECT_EQ( ids[id_pt_1] , id_0 ) ;
  }
  // Second set
  const uint32_t id_1 = ids[ nb_pts_1 ] ;
  for( int id_pt_2 = nb_pts_1 ; id_pt_2 < nb_pts_1 + nb_pts_2 ; ++id_pt_2 )
  {
    EXPECT_EQ( ids[id_pt_2] , id_1 ) ;
  }
  // Third set
  const uint32_t id_2 = ids[ nb_pts_1 + nb_pts_2 ] ;
  for( int id_pt_3 = nb_pts_1 + nb_pts_2 ; id_pt_3 < ids.size() ; ++id_pt_3 )
  {
    EXPECT_EQ( ids[id_pt_3] , id_2 ) ;
  }
}

TEST( clustering , threeClustersVec3 )
{
  std::vector< Vec3 > pts ;

  const Vec3 center1( -5.0 , -5.0 , -5.0 ) ;
  const Vec3 center2( 5.0 , 5.0 , 5.0 ) ;
  const Vec3 center3( 0.0 , 0.0 , 0.0 ) ;

  const int nb_pts_1 = 10000 ;
  const int nb_pts_2 = 10000 ;
  const int nb_pts_3 = 10000 ;

  std::mt19937_64 rng ;
  std::uniform_real_distribution<double> distrib( -2.0 , 2.0 ) ;

  // Generate points on each centers
  for( int id_pt_1 = 0 ; id_pt_1 < nb_pts_1 ; ++id_pt_1 )
  {
    const double dx = distrib( rng ) ;
    const double dy = distrib( rng ) ;
    const double dz = distrib( rng ) ;

    pts.push_back( center1 + Vec3( dx , dy , dz ) ) ;
  }
  for( int id_pt_2 = 0 ; id_pt_2 < nb_pts_2 ; ++id_pt_2 )
  {
    const double dx = distrib( rng ) ;
    const double dy = distrib( rng ) ;
    const double dz = distrib( rng ) ;

    pts.push_back( center2 + Vec3( dx , dy , dz  ) ) ;
  }
  for( int id_pt_3 = 0 ; id_pt_3 < nb_pts_3 ; ++id_pt_3 )
  {
    const double dx = distrib( rng ) ;
    const double dy = distrib( rng ) ;
    const double dz = distrib( rng ) ;

    pts.push_back( center3 + Vec3( dx , dy , dz ) ) ;
  }

  std::vector< uint32_t > ids ;
  std::vector< Vec3 > centers ;
  KMeans( pts , ids , centers , 3 ) ;

  std::cout << "Centers : " << std::endl ;
  for( const auto & it : centers )
  {
    std::cout << it << std::endl ;
  }

  // First set
  const uint32_t id_0 = ids[ 0 ] ;
  for( int id_pt_1 = 1 ; id_pt_1 < nb_pts_1 ; ++id_pt_1 )
  {
    EXPECT_EQ( ids[id_pt_1] , id_0 ) ;
  }
  // Second set
  const uint32_t id_1 = ids[ nb_pts_1 ] ;
  for( int id_pt_2 = nb_pts_1 ; id_pt_2 < nb_pts_1 + nb_pts_2 ; ++id_pt_2 )
  {
    EXPECT_EQ( ids[id_pt_2] , id_1 ) ;
  }
  // Third set
  const uint32_t id_2 = ids[ nb_pts_1 + nb_pts_2 ] ;
  for( int id_pt_3 = nb_pts_1 + nb_pts_2 ; id_pt_3 < ids.size() ; ++id_pt_3 )
  {
    EXPECT_EQ( ids[id_pt_3] , id_2 ) ;
  }
}


TEST( clustering , threeClustersStdArray )
{
  std::vector< std::array<double, 5> > pts ;

  const std::array<double, 5> center1 = { -5.0 , -5.0 , -5.0 , -5.0 , -5.0 } ;
  const std::array<double, 5> center2 = { 5.0 , 5.0 , 5.0 , 5.0 , 5.0 } ;
  const std::array<double, 5> center3 = { 0.0 , 0.0 , 0.0 , 0.0 , 0.0 } ;

  const int nb_pts_1 = 10000 ;
  const int nb_pts_2 = 10000 ;
  const int nb_pts_3 = 10000 ;

  std::mt19937_64 rng ;
  std::uniform_real_distribution<double> distrib( -2.0 , 2.0 ) ;

  // Generate points on each centers
  for( int id_pt_1 = 0 ; id_pt_1 < nb_pts_1 ; ++id_pt_1 )
  {
    const double dx = distrib( rng ) ;
    const double dy = distrib( rng ) ;
    const double dz = distrib( rng ) ;
    const double dw = distrib( rng ) ;
    const double dt = distrib( rng ) ;

    std::array<double, 5> tmp = center1 ;
    tmp[0] += dx ;
    tmp[1] += dy ;
    tmp[2] += dz ;
    tmp[3] += dw ;
    tmp[4] += dt ;
    pts.push_back( tmp ) ;
  }
  for( int id_pt_2 = 0 ; id_pt_2 < nb_pts_2 ; ++id_pt_2 )
  {
    const double dx = distrib( rng ) ;
    const double dy = distrib( rng ) ;
    const double dz = distrib( rng ) ;
    const double dw = distrib( rng ) ;
    const double dt = distrib( rng ) ;

    std::array<double, 5> tmp = center2 ;
    tmp[0] += dx ;
    tmp[1] += dy ;
    tmp[2] += dz ;
    tmp[3] += dw ;
    tmp[4] += dt ;
    pts.push_back( tmp ) ;
  }
  for( int id_pt_3 = 0 ; id_pt_3 < nb_pts_3 ; ++id_pt_3 )
  {
    const double dx = distrib( rng ) ;
    const double dy = distrib( rng ) ;
    const double dz = distrib( rng ) ;
    const double dw = distrib( rng ) ;
    const double dt = distrib( rng ) ;

    std::array<double, 5> tmp = center3 ;
    tmp[0] += dx ;
    tmp[1] += dy ;
    tmp[2] += dz ;
    tmp[3] += dw ;
    tmp[4] += dt ;
    pts.push_back( tmp ) ;
  }

  std::vector< uint32_t > ids ;
  std::vector< std::array<double, 5> > centers ;
  KMeans( pts , ids , centers , 3 ) ;

  std::cout << "Centers : " << std::endl ;
  for( const auto & it : centers )
  {
    for( const auto & val : it )
    {
      std::cout << val << " " ;
    }
    std::cout << std::endl ;
  }

  // First set
  const uint32_t id_0 = ids[ 0 ] ;
  for( int id_pt_1 = 1 ; id_pt_1 < nb_pts_1 ; ++id_pt_1 )
  {
    EXPECT_EQ( ids[id_pt_1] , id_0 ) ;
  }
  // Second set
  const uint32_t id_1 = ids[ nb_pts_1 ] ;
  for( int id_pt_2 = nb_pts_1 ; id_pt_2 < nb_pts_1 + nb_pts_2 ; ++id_pt_2 )
  {
    EXPECT_EQ( ids[id_pt_2] , id_1 ) ;
  }
  // Third set
  const uint32_t id_2 = ids[ nb_pts_1 + nb_pts_2 ] ;
  for( int id_pt_3 = nb_pts_1 + nb_pts_2 ; id_pt_3 < ids.size() ; ++id_pt_3 )
  {
    EXPECT_EQ( ids[id_pt_3] , id_2 ) ;
  }
}

TEST( clustering , threeClustersStdVector )
{
  std::vector< std::vector<double> > pts ;

  const std::vector<double> center1 = { -5.0 , -5.0 , -5.0 , -5.0 , -5.0 , -5.0 } ;
  const std::vector<double> center2 = { 5.0 , 5.0 , 5.0 , 5.0 , 5.0 , 5.0 } ;
  const std::vector<double> center3 = { 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 } ;

  const int nb_pts_1 = 10000 ;
  const int nb_pts_2 = 10000 ;
  const int nb_pts_3 = 10000 ;

  std::mt19937_64 rng ;
  std::uniform_real_distribution<double> distrib( -2.0 , 2.0 ) ;

  // Generate points on each centers
  for( int id_pt_1 = 0 ; id_pt_1 < nb_pts_1 ; ++id_pt_1 )
  {
    const double dx = distrib( rng ) ;
    const double dy = distrib( rng ) ;
    const double dz = distrib( rng ) ;
    const double dw = distrib( rng ) ;
    const double dt = distrib( rng ) ;
    const double ds = distrib( rng ) ;

    std::vector<double> tmp = center1 ;
    tmp[0] += dx ;
    tmp[1] += dy ;
    tmp[2] += dz ;
    tmp[3] += dw ;
    tmp[4] += dt ;
    tmp[5] += ds ;

    pts.push_back( tmp ) ;
  }
  for( int id_pt_2 = 0 ; id_pt_2 < nb_pts_2 ; ++id_pt_2 )
  {
    const double dx = distrib( rng ) ;
    const double dy = distrib( rng ) ;
    const double dz = distrib( rng ) ;
    const double dw = distrib( rng ) ;
    const double dt = distrib( rng ) ;
    const double ds = distrib( rng ) ;

    std::vector<double> tmp = center2 ;
    tmp[0] += dx ;
    tmp[1] += dy ;
    tmp[2] += dz ;
    tmp[3] += dw ;
    tmp[4] += dt ;
    tmp[5] += ds ;

    pts.push_back( tmp ) ;
  }
  for( int id_pt_3 = 0 ; id_pt_3 < nb_pts_3 ; ++id_pt_3 )
  {
    const double dx = distrib( rng ) ;
    const double dy = distrib( rng ) ;
    const double dz = distrib( rng ) ;
    const double dw = distrib( rng ) ;
    const double dt = distrib( rng ) ;
    const double ds = distrib( rng ) ;

    std::vector<double> tmp = center3 ;
    tmp[0] += dx ;
    tmp[1] += dy ;
    tmp[2] += dz ;
    tmp[3] += dw ;
    tmp[4] += dt ;
    tmp[5] += ds ;

    pts.push_back( tmp ) ;
  }

  std::vector< uint32_t > ids ;
  std::vector< std::vector<double> > centers ;
  KMeans( pts , ids , centers , 3 ) ;

  std::cout << "Centers : " << std::endl ;
  for( const auto & it : centers )
  {
    for( const auto & val : it )
    {
      std::cout << val << " " ;
    }
    std::cout << std::endl ;
  }

  // First set
  const uint32_t id_0 = ids[ 0 ] ;
  for( int id_pt_1 = 1 ; id_pt_1 < nb_pts_1 ; ++id_pt_1 )
  {
    EXPECT_EQ( ids[id_pt_1] , id_0 ) ;
  }
  // Second set
  const uint32_t id_1 = ids[ nb_pts_1 ] ;
  for( int id_pt_2 = nb_pts_1 ; id_pt_2 < nb_pts_1 + nb_pts_2 ; ++id_pt_2 )
  {
    EXPECT_EQ( ids[id_pt_2] , id_1 ) ;
  }
  // Third set
  const uint32_t id_2 = ids[ nb_pts_1 + nb_pts_2 ] ;
  for( int id_pt_3 = nb_pts_1 + nb_pts_2 ; id_pt_3 < ids.size() ; ++id_pt_3 )
  {
    EXPECT_EQ( ids[id_pt_3] , id_2 ) ;
  }
}


TEST( clustering , threeClustersEigenVec )
{
  std::vector< Vec > pts ;

  Vec center1( 6 ) , center2( 6 ) , center3( 6 ) ;
  center1 << -5.0 , -5.0 , -5.0 , -5.0 , -5.0 , -5.0  ;
  center2 <<  5.0 , 5.0 , 5.0 , 5.0 , 5.0 , 5.0 ;
  center3 <<  0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 ;

  const int nb_pts_1 = 10000 ;
  const int nb_pts_2 = 10000 ;
  const int nb_pts_3 = 10000 ;

  std::mt19937_64 rng ;
  std::uniform_real_distribution<double> distrib( -2.0 , 2.0 ) ;

  // Generate points on each centers
  for( int id_pt_1 = 0 ; id_pt_1 < nb_pts_1 ; ++id_pt_1 )
  {
    const double dx = distrib( rng ) ;
    const double dy = distrib( rng ) ;
    const double dz = distrib( rng ) ;
    const double dw = distrib( rng ) ;
    const double dt = distrib( rng ) ;
    const double ds = distrib( rng ) ;

    Vec tmp = center1 ;
    tmp[0] += dx ;
    tmp[1] += dy ;
    tmp[2] += dz ;
    tmp[3] += dw ;
    tmp[4] += dt ;
    tmp[5] += ds ;

    pts.push_back( tmp ) ;
  }
  for( int id_pt_2 = 0 ; id_pt_2 < nb_pts_2 ; ++id_pt_2 )
  {
    const double dx = distrib( rng ) ;
    const double dy = distrib( rng ) ;
    const double dz = distrib( rng ) ;
    const double dw = distrib( rng ) ;
    const double dt = distrib( rng ) ;
    const double ds = distrib( rng ) ;

    Vec tmp = center2 ;
    tmp[0] += dx ;
    tmp[1] += dy ;
    tmp[2] += dz ;
    tmp[3] += dw ;
    tmp[4] += dt ;
    tmp[5] += ds ;

    pts.push_back( tmp ) ;
  }
  for( int id_pt_3 = 0 ; id_pt_3 < nb_pts_3 ; ++id_pt_3 )
  {
    const double dx = distrib( rng ) ;
    const double dy = distrib( rng ) ;
    const double dz = distrib( rng ) ;
    const double dw = distrib( rng ) ;
    const double dt = distrib( rng ) ;
    const double ds = distrib( rng ) ;

    Vec tmp = center3 ;
    tmp[0] += dx ;
    tmp[1] += dy ;
    tmp[2] += dz ;
    tmp[3] += dw ;
    tmp[4] += dt ;
    tmp[5] += ds ;

    pts.push_back( tmp ) ;
  }

  std::vector< uint32_t > ids ;
  std::vector< Vec > centers ;
  KMeans( pts , ids , centers , 3 ) ;

  std::cout << "Centers : " << std::endl ;
  for( const auto & it : centers )
  {
    std::cout << it << std::endl ;
  }

  // First set
  const uint32_t id_0 = ids[ 0 ] ;
  for( int id_pt_1 = 1 ; id_pt_1 < nb_pts_1 ; ++id_pt_1 )
  {
    EXPECT_EQ( ids[id_pt_1] , id_0 ) ;
  }
  // Second set
  const uint32_t id_1 = ids[ nb_pts_1 ] ;
  for( int id_pt_2 = nb_pts_1 ; id_pt_2 < nb_pts_1 + nb_pts_2 ; ++id_pt_2 )
  {
    EXPECT_EQ( ids[id_pt_2] , id_1 ) ;
  }
  // Third set
  const uint32_t id_2 = ids[ nb_pts_1 + nb_pts_2 ] ;
  for( int id_pt_3 = nb_pts_1 + nb_pts_2 ; id_pt_3 < ids.size() ; ++id_pt_3 )
  {
    EXPECT_EQ( ids[id_pt_3] , id_2 ) ;
  }
}


/* ************************************************************************* */
int main()
{
  TestResult tr;
  return TestRegistry::runAllTests( tr );
}
/* ************************************************************************* */

