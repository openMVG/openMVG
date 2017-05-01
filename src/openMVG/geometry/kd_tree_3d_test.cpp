// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Romuald Perrot.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/geometry/kd_tree_3d.hpp"

#include "CppUnitLite/TestHarness.h"
#include "testing/testing.h"

#include <iostream>

using namespace openMVG;
using namespace openMVG::geometry;

TEST( kd_tree_3d, nothing )
{
  Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> data;

  KDTree3d<double> tree( data );

  Vec3 pt;
  int index;
  KDTree3d<double>::DistanceType dist;

  EXPECT_FALSE( tree.Search( pt, index, dist ) );
}

TEST( Kd_tree_3d, onePointOnePoint )
{
  Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> data;

  data.resize( 1, 3 );
  data( 0, 0 ) = 0.0;
  data( 0, 1 ) = 1.0;
  data( 0, 2 ) = 2.0;

  KDTree3d<double> tree( data );

  Vec3 pt( 0.0, 0.0, 0.0 );
  int index;
  KDTree3d<double>::DistanceType dist;

  EXPECT_TRUE( tree.Search( pt , index , dist ) ) ;

  EXPECT_EQ( index , 0 ) ;

  // squared distance 
  EXPECT_EQ( dist , 5.0 ) ;
}

TEST( Kd_tree_3d , ThreePointsOnePoint )
{
  Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> data;

  data.resize( 3, 3 );
  data( 0, 0 ) = 0.0;
  data( 0, 1 ) = 1.0;
  data( 0, 2 ) = 2.0;

  data( 1 , 0 ) = -10.0 ;
  data( 1 , 1 ) = -11.0 ;
  data( 1 , 2 ) = -3.0 ;

  data( 2 , 0 ) = 2.0 ;
  data( 2 , 1 ) = 3.0 ;
  data( 2 , 2 ) = 3.0 ; 

  KDTree3d<double> tree( data );

  Vec3 pt( 0.0, 0.0, 0.0 );
  int index;
  KDTree3d<double>::DistanceType dist;

  EXPECT_TRUE( tree.Search( pt , index , dist ) ) ;

  EXPECT_EQ( index , 0 ) ;

  // squared distance 
  EXPECT_EQ( dist , 5.0 ) ;


  pt = Vec3( 3.0 , 3.0 , 3.0 ) ;

  EXPECT_TRUE( tree.Search( pt , index , dist ) ) ;

  EXPECT_EQ( index , 2 ) ;

  // squared distance 
  EXPECT_EQ( dist , 1.0 ) ;


  pt = Vec3( -10.0 , -11.0 , 0.0 ) ;


  EXPECT_TRUE( tree.Search( pt , index , dist ) ) ;

  EXPECT_EQ( index , 1 ) ;

  // squared distance 
  EXPECT_EQ( dist , 9.0 ) ;
}


TEST( Kd_tree_3d , ThreePointsOnePointTwoNN )
{
  Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> data;

  data.resize( 3, 3 );
  data( 0, 0 ) = 0.0;
  data( 0, 1 ) = 1.0;
  data( 0, 2 ) = 2.0;

  data( 1 , 0 ) = -10.0 ;
  data( 1 , 1 ) = -11.0 ;
  data( 1 , 2 ) = -3.0 ;

  data( 2 , 0 ) = 2.0 ;
  data( 2 , 1 ) = 3.0 ;
  data( 2 , 2 ) = 3.0 ; 

  KDTree3d<double> tree( data );

  Vec3 pt( 0.0, 0.0, 0.0 );
  std::vector<int> index;
  std::vector<KDTree3d<double>::DistanceType> dist;

  index.resize(2) ;
  dist.resize(2) ;

  EXPECT_TRUE( tree.Search( pt , 2 , index , dist ) ) ;

  EXPECT_EQ( index[0] , 0 ) ;

  // squared distance 
  EXPECT_EQ( dist[0] , 5.0 ) ;

  EXPECT_EQ( index[1] , 2 ) ;

  EXPECT_EQ( dist[1] , 22.0 ) ;
}

TEST( Kd_tree_3d , ThreePointsTwoPoints )
{
  Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> data;

  data.resize( 3, 3 );
  data( 0, 0 ) = 0.0;
  data( 0, 1 ) = 1.0;
  data( 0, 2 ) = 2.0;

  data( 1 , 0 ) = -10.0 ;
  data( 1 , 1 ) = -11.0 ;
  data( 1 , 2 ) = -3.0 ;

  data( 2 , 0 ) = 2.0 ;
  data( 2 , 1 ) = 3.0 ;
  data( 2 , 2 ) = 3.0 ; 

  KDTree3d<double> tree( data );

  Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> pts ;

  pts.resize( 2 , 3 ) ;

  pts( 0 , 0 ) = 0.0 ;
  pts( 0 , 1 ) = 0.0 ; 
  pts( 0 , 2 ) = 0.0 ; 

  pts( 1 , 0 ) = -10.0 ;
  pts( 1 , 1 ) = -11.0 ;
  pts( 1 , 2 ) = 0.0 ; 

//  Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> indices ; 
//  Eigen::Matrix<KDTree3d<double>::DistanceType,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> dists ; 

  std::vector< int > indices ; 
  std::vector< KDTree3d<double>::DistanceType > dists ; 

  EXPECT_TRUE( tree.Search( pts , indices , dists ) ) ;

  // First pt 
  EXPECT_EQ( indices[0] , 0 ) ;
  EXPECT_EQ( dists[0] , 5.0 ) ;

  // Second pt 
  EXPECT_EQ( indices[1] , 1 ) ;
  EXPECT_EQ( dists[1] , 9.0 ) ;  
}

TEST( Kd_tree_3d , ThreePointsTwoPointsTwoQueries )
{
  Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> data;

  data.resize( 3, 3 );
  data( 0, 0 ) = 0.0;
  data( 0, 1 ) = 1.0;
  data( 0, 2 ) = 2.0;

  data( 1 , 0 ) = -10.0 ;
  data( 1 , 1 ) = -11.0 ;
  data( 1 , 2 ) = -3.0 ;

  data( 2 , 0 ) = 2.0 ;
  data( 2 , 1 ) = 3.0 ;
  data( 2 , 2 ) = 3.0 ; 

  KDTree3d<double> tree( data );

  Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> pts ;

  pts.resize( 2 , 3 ) ;

  pts( 0 , 0 ) = 0.0 ;
  pts( 0 , 1 ) = 0.0 ; 
  pts( 0 , 2 ) = 0.0 ; 

  pts( 1 , 0 ) = -10.0 ;
  pts( 1 , 1 ) = -11.0 ;
  pts( 1 , 2 ) = 0.0 ; 

  Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> indices ; 
  Eigen::Matrix<KDTree3d<double>::DistanceType,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> dists ; 

  EXPECT_TRUE( tree.Search( pts , 2 , indices , dists ) ) ;

  // First pt 
  EXPECT_EQ( indices(0,0) , 0 ) ;
  EXPECT_EQ( dists(0,0) , 5.0 ) ;
  EXPECT_EQ( indices(0,1) , 2 ) ;
  EXPECT_EQ( dists(0,1) , 22.0 ) ; 

  // Second pt 
  EXPECT_EQ( indices(1,0) , 1 ) ;
  EXPECT_EQ( dists(1,0) , 9.0 ) ;  
  EXPECT_EQ( indices(1,1) , 0 ) ;
  EXPECT_EQ( dists(1,1) , 248.0 )
  
}


/* ************************************************************************* */
int main()
{
  TestResult tr;
  return TestRegistry::runAllTests( tr );
}
/* ************************************************************************* */
