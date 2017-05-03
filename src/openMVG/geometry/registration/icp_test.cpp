// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Romuald Perrot.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "testing/testing.h"

#include "openMVG/geometry/registration/icp.hpp"

#include <cmath>

TEST( icp , icp_no_transfo_f )
{
  // Generate point set
  Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor> target;
  Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor> data;

  // A sphere
  const double pi          = 3.14159265358979;
  const int nb_x = 64 ;
  const int nb_y = 64 ;

  const double x_min = -10 ;
  const double x_max = 10 ;
  const double y_min = -10 ;
  const double y_max = 10 ;

  const double x_range = x_max - x_min ;
  const double y_range = y_max - y_min ;
  const double delta_x = x_range / static_cast<double>( nb_x ) ;
  const double delta_y = y_range / static_cast<double>( nb_y ) ;

  target.resize( nb_x * nb_y , 3 );
  data.resize( nb_x * nb_y , 3 );

  int id_pt = 0;
  for( int id_y = 0 ; id_y < nb_y ; ++id_y )
  {
    for( int id_x = 0 ; id_x < nb_x ; ++id_x )
    {
      const double x = x_min + id_x * delta_x ;
      const double y = y_min + id_y * delta_y ;

      const double rad = std::sqrt( x * x + y * y ) ;
      const double z = std::max( 0.0 , x ) + std::max( 0.0 , y ) + std::cos( 2.0 * rad ) ;

      target( id_pt, 0 ) = x;
      target( id_pt, 1 ) = y;
      target( id_pt, 2 ) = z;

      data( id_pt, 0 ) = x ;
      data( id_pt, 1 ) = y ;
      data( id_pt, 2 ) = z ;

      ++id_pt;
    }
  }

  openMVG::Mat3 R;
  openMVG::Vec3 t;

  openMVG::geometry::registration::ICP( target, data, 200, 0.0001f , t, R );
  const openMVG::Vec3 invT = R.transpose() * t ;

  std::cout << "Computed" << std::endl;
  std::cout << invT << std::endl;
  std::cout << R << std::endl;

  const double tx = 0.0 ;
  const double ty = 0.0 ;
  const double tz = 0.0 ;

  openMVG::Mat3 r ;
  r.setIdentity() ;
  std::cout << "Expected" << std::endl;
  std::cout << -tx << " - " << -ty << " - " << -tz << std::endl;
  std::cout << r << std::endl ;

  EXPECT_NEAR( invT[ 0 ], -tx, 0.00001 );
  EXPECT_NEAR( invT[ 1 ], -ty, 0.00001 );
  EXPECT_NEAR( invT[ 2 ], -tz, 0.00001 );

  EXPECT_NEAR( R( 0, 0 ), r( 0, 0 ), 0.00001 );
  EXPECT_NEAR( R( 0, 1 ), r( 0, 1 ), 0.00001 );
  EXPECT_NEAR( R( 0, 2 ), r( 0, 2 ), 0.00001 );

  EXPECT_NEAR( R( 1, 0 ), r( 1, 0 ), 0.00001 );
  EXPECT_NEAR( R( 1, 1 ), r( 1, 1 ), 0.00001 );
  EXPECT_NEAR( R( 1, 2 ), r( 1, 2 ), 0.00001 );

  EXPECT_NEAR( R( 2, 0 ), r( 2, 0 ), 0.00001 );
  EXPECT_NEAR( R( 2, 1 ), r( 2, 1 ), 0.00001 );
  EXPECT_NEAR( R( 2, 2 ), r( 2, 2 ), 0.00001 );
}

TEST( icp, icp_no_transfo )
{
  // Generate point set
  Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> target;
  Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> data;

  // A sphere
  const double pi          = 3.14159265358979;
  const int nb_x = 64 ;
  const int nb_y = 64 ;

  const double x_min = -10 ;
  const double x_max = 10 ;
  const double y_min = -10 ;
  const double y_max = 10 ;

  const double x_range = x_max - x_min ;
  const double y_range = y_max - y_min ;
  const double delta_x = x_range / static_cast<double>( nb_x ) ;
  const double delta_y = y_range / static_cast<double>( nb_y ) ;

  target.resize( nb_x * nb_y , 3 );
  data.resize( nb_x * nb_y , 3 );

  int id_pt = 0;
  for( int id_y = 0 ; id_y < nb_y ; ++id_y )
  {
    for( int id_x = 0 ; id_x < nb_x ; ++id_x )
    {
      const double x = x_min + id_x * delta_x ;
      const double y = y_min + id_y * delta_y ;

      const double rad = std::sqrt( x * x + y * y ) ;
      const double z = std::max( 0.0 , x ) + std::max( 0.0 , y ) + std::cos( 2.0 * rad ) ;

      target( id_pt, 0 ) = x;
      target( id_pt, 1 ) = y;
      target( id_pt, 2 ) = z;

      data( id_pt, 0 ) = x ;
      data( id_pt, 1 ) = y ;
      data( id_pt, 2 ) = z ;

      ++id_pt;
    }
  }

  openMVG::Mat3 R;
  openMVG::Vec3 t;

  openMVG::geometry::registration::ICP( target, data, 200, 0.0001, t, R );
  const openMVG::Vec3 invT = R.transpose() * t ;

  std::cout << "Computed" << std::endl;
  std::cout << invT << std::endl;
  std::cout << R << std::endl;

  const double tx = 0.0 ;
  const double ty = 0.0 ;
  const double tz = 0.0 ;

  openMVG::Mat3 r ;
  r.setIdentity() ;
  std::cout << "Expected" << std::endl;
  std::cout << -tx << " - " << -ty << " - " << -tz << std::endl;
  std::cout << r << std::endl ;

  EXPECT_NEAR( invT[ 0 ], -tx, 0.00001 );
  EXPECT_NEAR( invT[ 1 ], -ty, 0.00001 );
  EXPECT_NEAR( invT[ 2 ], -tz, 0.00001 );

  EXPECT_NEAR( R( 0, 0 ), r( 0, 0 ), 0.00001 );
  EXPECT_NEAR( R( 0, 1 ), r( 0, 1 ), 0.00001 );
  EXPECT_NEAR( R( 0, 2 ), r( 0, 2 ), 0.00001 );

  EXPECT_NEAR( R( 1, 0 ), r( 1, 0 ), 0.00001 );
  EXPECT_NEAR( R( 1, 1 ), r( 1, 1 ), 0.00001 );
  EXPECT_NEAR( R( 1, 2 ), r( 1, 2 ), 0.00001 );

  EXPECT_NEAR( R( 2, 0 ), r( 2, 0 ), 0.00001 );
  EXPECT_NEAR( R( 2, 1 ), r( 2, 1 ), 0.00001 );
  EXPECT_NEAR( R( 2, 2 ), r( 2, 2 ), 0.00001 );
}

TEST( icp, icp_rotate )
{
  // Generate point set
  Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> target;
  Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> data;

  // A sphere
  const double pi          = 3.14159265358979;
  openMVG::Mat3 r ;
  r = Eigen::AngleAxis<double>( pi / 10.0, openMVG::Vec3( 1, 1, 1 ).normalized() );

  const int nb_x = 64 ;
  const int nb_y = 64 ;

  const double x_min = -10 ;
  const double x_max = 10 ;
  const double y_min = -10 ;
  const double y_max = 10 ;

  const double x_range = x_max - x_min ;
  const double y_range = y_max - y_min ;
  const double delta_x = x_range / static_cast<double>( nb_x ) ;
  const double delta_y = y_range / static_cast<double>( nb_y ) ;

  target.resize( nb_x * nb_y , 3 );
  data.resize( nb_x * nb_y , 3 );

  int id_pt = 0;
  for( int id_y = 0 ; id_y < nb_y ; ++id_y )
  {
    for( int id_x = 0 ; id_x < nb_x ; ++id_x )
    {
      const double x = x_min + id_x * delta_x ;
      const double y = y_min + id_y * delta_y ;

      const double rad = std::sqrt( x * x + y * y ) ;
      const double z = std::max( 0.0 , x ) + std::max( 0.0 , y ) + std::cos( 2.0 * rad ) ;

      target( id_pt, 0 ) = x;
      target( id_pt, 1 ) = y;
      target( id_pt, 2 ) = z;

      const openMVG::Vec3 orig( x, y, z );
      const openMVG::Vec3 tra = r * orig;

      data( id_pt, 0 ) = tra[ 0 ];
      data( id_pt, 1 ) = tra[ 1 ];
      data( id_pt, 2 ) = tra[ 2 ];

      ++id_pt;
    }
  }

  // Inverse here
  r = r.inverse().eval();

  openMVG::Mat3 R;
  openMVG::Vec3 t;

  openMVG::geometry::registration::ICP( target, data, 200, 0.0001, t, R );
  const openMVG::Vec3 invT = R.transpose() * t ;

  std::cout << "Computed" << std::endl;
  std::cout << invT << std::endl;
  std::cout << R << std::endl;

  std::cout << "Expected" << std::endl;
  std::cout << r << std::endl;


  EXPECT_NEAR( t[ 0 ], 0, 0.00001 );
  EXPECT_NEAR( t[ 1 ], 0, 0.00001 );
  EXPECT_NEAR( t[ 2 ], 0, 0.00001 );

  EXPECT_NEAR( R( 0, 0 ), r( 0, 0 ), 0.00001 );
  EXPECT_NEAR( R( 0, 1 ), r( 0, 1 ), 0.00001 );
  EXPECT_NEAR( R( 0, 2 ), r( 0, 2 ), 0.00001 );

  EXPECT_NEAR( R( 1, 0 ), r( 1, 0 ), 0.00001 );
  EXPECT_NEAR( R( 1, 1 ), r( 1, 1 ), 0.00001 );
  EXPECT_NEAR( R( 1, 2 ), r( 1, 2 ), 0.00001 );

  EXPECT_NEAR( R( 2, 0 ), r( 2, 0 ), 0.00001 );
  EXPECT_NEAR( R( 2, 1 ), r( 2, 1 ), 0.00001 );
  EXPECT_NEAR( R( 2, 2 ), r( 2, 2 ), 0.00001 );
}

TEST( icp, icp_translate )
{
  // Generate point set
  Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> target;
  Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> data;

  const double tx = 10;
  const double ty = 12;
  const double tz = 10;

  const double pi          = 3.14159265358979;


  const int nb_x = 64 ;
  const int nb_y = 64 ;

  const double x_min = -10 ;
  const double x_max = 10 ;
  const double y_min = -10 ;
  const double y_max = 10 ;

  const double x_range = x_max - x_min ;
  const double y_range = y_max - y_min ;
  const double delta_x = x_range / static_cast<double>( nb_x ) ;
  const double delta_y = y_range / static_cast<double>( nb_y ) ;

  target.resize( nb_x * nb_y , 3 );
  data.resize( nb_x * nb_y , 3 );

  int id_pt = 0;
  for( int id_y = 0 ; id_y < nb_y ; ++id_y )
  {
    for( int id_x = 0 ; id_x < nb_x ; ++id_x )
    {
      const double x = x_min + id_x * delta_x ;
      const double y = y_min + id_y * delta_y ;

      const double rad = std::sqrt( x * x + y * y ) ;
      const double z = std::max( 0.0 , x ) + std::max( 0.0 , y ) + std::cos( 2.0 * rad ) ;

      target( id_pt, 0 ) = x;
      target( id_pt, 1 ) = y;
      target( id_pt, 2 ) = z;

      data( id_pt, 0 ) = x + tx;
      data( id_pt, 1 ) = y + ty;
      data( id_pt, 2 ) = z + tz;

      ++id_pt;
    }
  }

  openMVG::Mat3 R;
  openMVG::Vec3 t;

  openMVG::geometry::registration::ICP( target, data, 2000, 0.0001, t, R );
  const openMVG::Vec3 invT = R.transpose() * t ;

  std::cout << "Computed" << std::endl;
  std::cout << invT << std::endl;
  std::cout << R << std::endl;

  openMVG::Mat3 r ;
  r.setIdentity() ;
  std::cout << "Expected" << std::endl;
  std::cout << -tx << " - " << -ty << " - " << -tz << std::endl;
  std::cout << r << std::endl ;

  EXPECT_NEAR( invT[ 0 ], -tx, 0.00001 );
  EXPECT_NEAR( invT[ 1 ], -ty, 0.00001 );
  EXPECT_NEAR( invT[ 2 ], -tz, 0.00001 );

  EXPECT_NEAR( R( 0, 0 ), r( 0, 0 ), 0.00001 );
  EXPECT_NEAR( R( 0, 1 ), r( 0, 1 ), 0.00001 );
  EXPECT_NEAR( R( 0, 2 ), r( 0, 2 ), 0.00001 );

  EXPECT_NEAR( R( 1, 0 ), r( 1, 0 ), 0.00001 );
  EXPECT_NEAR( R( 1, 1 ), r( 1, 1 ), 0.00001 );
  EXPECT_NEAR( R( 1, 2 ), r( 1, 2 ), 0.00001 );

  EXPECT_NEAR( R( 2, 0 ), r( 2, 0 ), 0.00001 );
  EXPECT_NEAR( R( 2, 1 ), r( 2, 1 ), 0.00001 );
  EXPECT_NEAR( R( 2, 2 ), r( 2, 2 ), 0.00001 );
}

TEST( icp, icp_translate_normal )
{
  // Generate point set
  Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> target;
  Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> target_n;
  Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> data;

  const double tx = 10;
  const double ty = 12;
  const double tz = 10;

  const double pi          = 3.14159265358979;

  const int nb_x = 64 ;
  const int nb_y = 64 ;

  const double x_min = -10 ;
  const double x_max = 10 ;
  const double y_min = -10 ;
  const double y_max = 10 ;

  const double x_range = x_max - x_min ;
  const double y_range = y_max - y_min ;
  const double delta_x = x_range / static_cast<double>( nb_x ) ;
  const double delta_y = y_range / static_cast<double>( nb_y ) ;

  target.resize( nb_x * nb_y , 3 );
  target_n.resize( nb_x * nb_y , 3 ) ;
  data.resize( nb_x * nb_y , 3 );

  int id_pt = 0;
  for( int id_y = 0 ; id_y < nb_y ; ++id_y )
  {
    for( int id_x = 0 ; id_x < nb_x ; ++id_x )
    {
      const double x = x_min + id_x * delta_x ;
      const double y = y_min + id_y * delta_y ;

      const double rad = std::sqrt( x * x + y * y ) ;
      const double z = std::max( 0.0 , x ) + std::max( 0.0 , y ) + std::cos( 2.0 * rad ) ;

      const double dx_x = 1.0 ;
      const double dy_x = 0.0 ;
      const double dz_x = std::max( 0.0 , x ) / x - 2.0 * std::sin( 2.0 * rad ) / rad ;

      const double dx_y = 0.0 ;
      const double dy_y = 1.0 ;
      const double dz_y = std::max( 0.0 , y ) / y - 2.0 * std::sin( 2.0 * rad ) / rad ;

      const openMVG::Vec3 vx( dx_x , dx_y , dz_x ) ;
      const openMVG::Vec3 vy( dy_x , dy_y , dz_y ) ;

      const openMVG::Vec3 n = ( vx.cross( vy ) ).normalized() ;

      target( id_pt, 0 ) = x;
      target( id_pt, 1 ) = y;
      target( id_pt, 2 ) = z;

      target_n( id_pt , 0 ) = n[0] ;
      target_n( id_pt , 1 ) = n[1] ;
      target_n( id_pt , 2 ) = n[2] ;

      data( id_pt, 0 ) = x + tx;
      data( id_pt, 1 ) = y + ty;
      data( id_pt, 2 ) = z + tz;

      ++id_pt;
    }
  }

  openMVG::Mat3 R;
  openMVG::Vec3 t;

  openMVG::geometry::registration::ICP( target, data, 2000, 0.0001, t, R );
  const openMVG::Vec3 invT = R.transpose() * t ;

  std::cout << "Computed" << std::endl;
  std::cout << invT << std::endl;
  std::cout << R << std::endl;

  openMVG::Mat3 r ;
  r.setIdentity() ;
  std::cout << "Expected" << std::endl;
  std::cout << -tx << " - " << -ty << " - " << -tz << std::endl;
  std::cout << r << std::endl ;


  EXPECT_NEAR( invT[ 0 ], -tx, 0.00001 );
  EXPECT_NEAR( invT[ 1 ], -ty, 0.00001 );
  EXPECT_NEAR( invT[ 2 ], -tz, 0.00001 );

  EXPECT_NEAR( R( 0, 0 ), r( 0, 0 ), 0.00001 );
  EXPECT_NEAR( R( 0, 1 ), r( 0, 1 ), 0.00001 );
  EXPECT_NEAR( R( 0, 2 ), r( 0, 2 ), 0.00001 );

  EXPECT_NEAR( R( 1, 0 ), r( 1, 0 ), 0.00001 );
  EXPECT_NEAR( R( 1, 1 ), r( 1, 1 ), 0.00001 );
  EXPECT_NEAR( R( 1, 2 ), r( 1, 2 ), 0.00001 );

  EXPECT_NEAR( R( 2, 0 ), r( 2, 0 ), 0.00001 );
  EXPECT_NEAR( R( 2, 1 ), r( 2, 1 ), 0.00001 );
  EXPECT_NEAR( R( 2, 2 ), r( 2, 2 ), 0.00001 );
}

TEST( icp, icp_rot_translate )
{
  // Generate point set
  Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> target;
  Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> data;

  const double pi          = 3.14159265358979;
  openMVG::Mat3 r ;
  r = Eigen::AngleAxis<double>( pi / 10.0, openMVG::Vec3( 1, 1, 1 ).normalized() );

  const double tx = 0.1;
  const double ty = 0.2;
  const double tz = 0.3;

  const int nb_x = 64 ;
  const int nb_y = 64 ;

  const double x_min = -10 ;
  const double x_max = 10 ;
  const double y_min = -10 ;
  const double y_max = 10 ;

  const double x_range = x_max - x_min ;
  const double y_range = y_max - y_min ;
  const double delta_x = x_range / static_cast<double>( nb_x ) ;
  const double delta_y = y_range / static_cast<double>( nb_y ) ;

  target.resize( nb_x * nb_y , 3 );
  data.resize( nb_x * nb_y , 3 );

  int id_pt = 0;
  for( int id_y = 0 ; id_y < nb_y ; ++id_y )
  {
    for( int id_x = 0 ; id_x < nb_x ; ++id_x )
    {
      const double x = x_min + id_x * delta_x ;
      const double y = y_min + id_y * delta_y ;

      const double rad = std::sqrt( x * x + y * y ) ;
      const double z = std::max( 0.0 , x ) + std::max( 0.0 , y ) + std::cos( 2.0 * rad ) ;

      target( id_pt, 0 ) = x;
      target( id_pt, 1 ) = y;
      target( id_pt, 2 ) = z;

      const openMVG::Vec3 orig( x, y, z );
      const openMVG::Vec3 tra = r * orig;

      data( id_pt, 0 ) = tra[0] + tx;
      data( id_pt, 1 ) = tra[1] + ty;
      data( id_pt, 2 ) = tra[2] + tz;

      ++id_pt;
    }
  }
  std::cout << Eigen::umeyama( data.transpose() , target.transpose() ) ;

  // Inverse here
  r = r.inverse().eval();

  std::cout << r * openMVG::Vec3( tx , ty , tz ) ;

  openMVG::Mat3 R;
  openMVG::Vec3 t;

  openMVG::geometry::registration::ICP( target, data, 2000 , 0.00001 , t, R );

  const openMVG::Vec3 invT = R.transpose() * t ;

  std::cout << "Computed" << std::endl;
  std::cout << invT << std::endl;
  std::cout << R << std::endl;


  std::cout << "Expected" << std::endl;
  std::cout << -tx << " - " << -ty << " - " << -tz << std::endl;
  std::cout << r << std::endl;

  EXPECT_NEAR( invT[ 0 ], -tx, 0.00001 );
  EXPECT_NEAR( invT[ 1 ], -ty, 0.00001 );
  EXPECT_NEAR( invT[ 2 ], -tz, 0.00001 );

  EXPECT_NEAR( R( 0, 0 ), r( 0, 0 ), 0.00001 );
  EXPECT_NEAR( R( 0, 1 ), r( 0, 1 ), 0.00001 );
  EXPECT_NEAR( R( 0, 2 ), r( 0, 2 ), 0.00001 );

  EXPECT_NEAR( R( 1, 0 ), r( 1, 0 ), 0.00001 );
  EXPECT_NEAR( R( 1, 1 ), r( 1, 1 ), 0.00001 );
  EXPECT_NEAR( R( 1, 2 ), r( 1, 2 ), 0.00001 );

  EXPECT_NEAR( R( 2, 0 ), r( 2, 0 ), 0.00001 );
  EXPECT_NEAR( R( 2, 1 ), r( 2, 1 ), 0.00001 );
  EXPECT_NEAR( R( 2, 2 ), r( 2, 2 ), 0.00001 );
}


TEST( icp, icp_rot_translate_normal )
{
  // Generate point set
  Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> target;
  Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> target_n;
  Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> data;

  openMVG::Mat3 r;
  const double pi          = 3.14159265358979;
  r = Eigen::AngleAxis<double>( pi / 10.0, openMVG::Vec3( 1, 1, 1 ).normalized() );

  const double tx = 0.1;
  const double ty = 0.2;
  const double tz = 0.3;

  const int nb_x = 64 ;
  const int nb_y = 64 ;

  const double x_min = -10 ;
  const double x_max = 10 ;
  const double y_min = -10 ;
  const double y_max = 10 ;

  const double x_range = x_max - x_min ;
  const double y_range = y_max - y_min ;
  const double delta_x = x_range / static_cast<double>( nb_x ) ;
  const double delta_y = y_range / static_cast<double>( nb_y ) ;

  target.resize( nb_x * nb_y , 3 );
  target_n.resize( nb_x * nb_y , 3 ) ;
  data.resize( nb_x * nb_y , 3 );

  int id_pt = 0;
  for( int id_y = 0 ; id_y < nb_y ; ++id_y )
  {
    for( int id_x = 0 ; id_x < nb_x ; ++id_x )
    {
      const double x = x_min + id_x * delta_x ;
      const double y = y_min + id_y * delta_y ;

      const double rad = std::sqrt( x * x + y * y ) ;
      const double z = std::max( 0.0 , x ) + std::max( 0.0 , y ) + std::cos( 2.0 * rad ) ;

      const double dx_x = 1.0 ;
      const double dy_x = 0.0 ;
      const double dz_x = std::max( 0.0 , x ) / x - 2.0 * std::sin( 2.0 * rad ) / rad ;

      const double dx_y = 0.0 ;
      const double dy_y = 1.0 ;
      const double dz_y = std::max( 0.0 , y ) / y - 2.0 * std::sin( 2.0 * rad ) / rad ;

      const openMVG::Vec3 vx( dx_x , dx_y , dz_x ) ;
      const openMVG::Vec3 vy( dy_x , dy_y , dz_y ) ;

      const openMVG::Vec3 n = ( vx.cross( vy ) ).normalized() ;

      target( id_pt, 0 ) = x;
      target( id_pt, 1 ) = y;
      target( id_pt, 2 ) = z;

      target_n( id_pt , 0 ) = n[0] ;
      target_n( id_pt , 1 ) = n[1] ;
      target_n( id_pt , 2 ) = n[2] ;

      const openMVG::Vec3 orig( x, y, z );
      const openMVG::Vec3 tra = r * orig;

      data( id_pt, 0 ) = tra[0] + tx;
      data( id_pt, 1 ) = tra[1] + ty;
      data( id_pt, 2 ) = tra[2] + tz;

      ++id_pt;
    }
  }

  // Inverse here
  r = r.inverse().eval();

  openMVG::Mat3 R;
  openMVG::Vec3 t;

  openMVG::geometry::registration::ICP( target, data, 2000 , 0.0001 , t, R );

  const openMVG::Vec3 invT = R.transpose() * t ;

  std::cout << "Computed" << std::endl;
  std::cout << invT << std::endl;
  std::cout << R << std::endl;

  std::cout << "Expected" << std::endl;
  std::cout << -tx << " - " << -ty << " - " << -tz << std::endl;
  std::cout << r << std::endl;

  EXPECT_NEAR( invT[ 0 ], -tx, 0.00001 );
  EXPECT_NEAR( invT[ 1 ], -ty, 0.00001 );
  EXPECT_NEAR( invT[ 2 ], -tz, 0.00001 );

  EXPECT_NEAR( R( 0, 0 ), r( 0, 0 ), 0.00001 );
  EXPECT_NEAR( R( 0, 1 ), r( 0, 1 ), 0.00001 );
  EXPECT_NEAR( R( 0, 2 ), r( 0, 2 ), 0.00001 );

  EXPECT_NEAR( R( 1, 0 ), r( 1, 0 ), 0.00001 );
  EXPECT_NEAR( R( 1, 1 ), r( 1, 1 ), 0.00001 );
  EXPECT_NEAR( R( 1, 2 ), r( 1, 2 ), 0.00001 );

  EXPECT_NEAR( R( 2, 0 ), r( 2, 0 ), 0.00001 );
  EXPECT_NEAR( R( 2, 1 ), r( 2, 1 ), 0.00001 );
  EXPECT_NEAR( R( 2, 2 ), r( 2, 2 ), 0.00001 );
}


/* ************************************************************************* */
int main()
{
  TestResult tr;
  return TestRegistry::runAllTests( tr );
}
/* ************************************************************************* */
