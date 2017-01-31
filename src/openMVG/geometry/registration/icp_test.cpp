#include "testing/testing.h"

#include "openMVG/geometry/registration/icp.hpp"

#include <cmath>

TEST( icp , icp_no_transfo )
{
  // Generate point set 
  Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> target ;
  Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> data ; 

  // A sphere 
  const double pi = 3.14159265358979 ;
  const int nb_theta = 32 ;
  const int nb_phi   = 16 ;
  const double delta_theta = ( 2.0 * pi ) / nb_theta ;
  const double delta_phi = ( pi ) / nb_phi ;
  const double rad = 2.0 ;

  int nb_pt = 0 ;
  for( double theta = 0.0 ; theta < 2.0 * pi ; theta += delta_theta )
  {
    for( double phi = -pi/2.0 ; phi < pi/2.0 ; phi += delta_phi )
    {
      ++nb_pt ; 
    }
  }
  
  target.resize( nb_pt , 3 ) ;
  data.resize( nb_pt , 3 ) ; 

  int id_pt = 0 ; 
  for( double theta = 0.0 ; theta < 2.0 * pi ; theta += delta_theta )
  {
    for( double phi = -pi/2.0 ; phi < pi/2.0 ; phi += delta_phi )
    {
      const double x = rad * std::cos( phi ) * std::cos( theta ) ;
      const double y = rad * std::cos( phi ) * std::sin( theta ) ;
      const double z = rad * std::sin( phi ) ;
    
      target( id_pt , 0 ) = x ;
      target( id_pt , 1 ) = y ;
      target( id_pt , 2 ) = z ; 

      data( id_pt , 0 ) = x ; 
      data( id_pt , 1 ) = y ;
      data( id_pt , 2 ) = z ; 

      ++id_pt ;   
    }
  }

  openMVG::Mat3 R ;
  openMVG::Vec3 t ; 

  openMVG::geometry::registration::ICP( target , data , 200 , 0.0 , t , R ) ; 

  EXPECT_EQ( t[0] , 0.0 ) ;
  EXPECT_EQ( t[1] , 0.0 ) ;
  EXPECT_EQ( t[2] , 0.0 ) ; 

  EXPECT_EQ( R(0,0) , 1.0 ) ;
  EXPECT_EQ( R(0,1) , 0.0 ) ;
  EXPECT_EQ( R(0,2) , 0.0 ) ;

  EXPECT_EQ( R(1,0) , 0.0 ) ;
  EXPECT_EQ( R(1,1) , 1.0 ) ;
  EXPECT_EQ( R(1,2) , 0.0 ) ;

  EXPECT_EQ( R(2,0) , 0.0 ) ;
  EXPECT_EQ( R(2,1) , 0.0 ) ;
  EXPECT_EQ( R(2,2) , 1.0 ) ; 
}


TEST( icp , icp_rotate )
{
  // Generate point set 
  Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> target ;
  Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> data ; 

  // A sphere 
  const double pi = 3.14159265358979 ;
  const int nb_theta = 32 ;
  const int nb_phi   = 16 ;
  const double delta_theta = ( 2.0 * pi ) / nb_theta ;
  const double delta_phi = ( pi ) / nb_phi ;
  const double rad = 2.0 ;

  int nb_pt = 0 ;
  for( double theta = 0.0 ; theta < 2.0 * pi ; theta += delta_theta )
  {
    for( double phi = -pi/2.0 ; phi < pi/2.0 ; phi += delta_phi )
    {
      ++nb_pt ; 
    }
  }
  
  target.resize( nb_pt , 3 ) ;
  data.resize( nb_pt , 3 ) ; 

  openMVG::Mat3 r ;
  r = Eigen::AngleAxis<double>( pi / 3.0 , openMVG::Vec3( 1 , 1 , 1 ).normalized() ) ;


  int id_pt = 0 ; 
  for( double theta = 0.0 ; theta < 2.0 * pi ; theta += delta_theta )
  {
    for( double phi = -pi/2.0 ; phi < pi/2.0 ; phi += delta_phi )
    {
      const double x = rad * std::cos( phi ) * std::cos( theta ) ;
      const double y = rad * std::cos( phi ) * std::sin( theta ) ;
      const double z = rad * std::sin( phi ) ;
    
      target( id_pt , 0 ) = x ;
      target( id_pt , 1 ) = y ;
      target( id_pt , 2 ) = z ; 


      const openMVG::Vec3 orig( x , y , z ) ;
      const openMVG::Vec3 tra = r * orig ; 
      
      data( id_pt , 0 ) = tra[0] ; 
      data( id_pt , 1 ) = tra[1] ;
      data( id_pt , 2 ) = tra[2] ; 

      ++id_pt ;   
    }
  }

  // Inverse here 
  r = r.inverse().eval() ; 

  openMVG::Mat3 R ;
  openMVG::Vec3 t ; 

  openMVG::geometry::registration::ICP( target , data , 1000 , 0.0 , t , R ) ; 

  EXPECT_NEAR( t[0] , 0 , 0.00001 ) ;
  EXPECT_NEAR( t[1] , 0 , 0.00001 ) ;
  EXPECT_NEAR( t[2] , 0 , 0.00001 ) ; 

  EXPECT_NEAR( R(0,0) , r(0,0) , 0.00001 ) ;
  EXPECT_NEAR( R(0,1) , r(0,1) , 0.00001 ) ;
  EXPECT_NEAR( R(0,2) , r(0,2) , 0.00001 ) ;

  EXPECT_NEAR( R(1,0) , r(1,0) , 0.00001 ) ;
  EXPECT_NEAR( R(1,1) , r(1,1) , 0.00001 ) ;
  EXPECT_NEAR( R(1,2) , r(1,2) , 0.00001 ) ;

  EXPECT_NEAR( R(2,0) , r(2,0) , 0.00001 ) ;
  EXPECT_NEAR( R(2,1) , r(2,1) , 0.00001 ) ;
  EXPECT_NEAR( R(2,2) , r(2,2) , 0.00001 ) ; 
}

TEST( icp , icp_translate )
{
  // Generate point set 
  Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> target ;
  Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> data ; 

  // A sphere 
  const double pi = 3.14159265358979 ;
  const int nb_theta = 32 ;
  const int nb_phi   = 16 ;
  const double delta_theta = ( 2.0 * pi ) / nb_theta ;
  const double delta_phi = ( pi ) / nb_phi ;
  const double rad = 2.0 ;

  int nb_pt = 0 ;
  for( double theta = 0.0 ; theta < 2.0 * pi ; theta += delta_theta )
  {
    for( double phi = -pi/2.0 ; phi < pi/2.0 ; phi += delta_phi )
    {
      ++nb_pt ; 
    }
  }
  
  target.resize( nb_pt , 3 ) ;
  data.resize( nb_pt , 3 ) ; 

  const double tx = 10.0 ;
  const double ty = 20.0 ;
  const double tz = 10.0 ; 

  int id_pt = 0 ; 
  for( double theta = 0.0 ; theta < 2.0 * pi ; theta += delta_theta )
  {
    for( double phi = -pi/2.0 ; phi < pi/2.0 ; phi += delta_phi )
    {
      const double x = rad * std::cos( phi ) * std::cos( theta ) ;
      const double y = rad * std::cos( phi ) * std::sin( theta ) ;
      const double z = rad * std::sin( phi ) ;
    
      target( id_pt , 0 ) = x ;
      target( id_pt , 1 ) = y ;
      target( id_pt , 2 ) = z ; 

      data( id_pt , 0 ) = x + tx ; 
      data( id_pt , 1 ) = y + ty ;
      data( id_pt , 2 ) = z + tz ; 

      ++id_pt ;   
    }
  }

  openMVG::Mat3 R ;
  openMVG::Vec3 t ; 

  openMVG::geometry::registration::ICP( target , data , 1000 , 0.0 , t , R ) ; 

  EXPECT_NEAR( t[0] , -tx , 0.00001 ) ;
  EXPECT_NEAR( t[1] , -ty , 0.00001 ) ;
  EXPECT_NEAR( t[2] , -tz , 0.00001 ) ; 

  EXPECT_NEAR( R(0,0) , 1.0 , 0.00001 ) ;
  EXPECT_NEAR( R(0,1) , 0.0 , 0.00001 ) ;
  EXPECT_NEAR( R(0,2) , 0.0 , 0.00001 ) ;

  EXPECT_NEAR( R(1,0) , 0.0 , 0.00001 ) ;
  EXPECT_NEAR( R(1,1) , 1.0 , 0.00001 ) ;
  EXPECT_NEAR( R(1,2) , 0.0 , 0.00001 ) ;

  EXPECT_NEAR( R(2,0) , 0.0 , 0.00001 ) ;
  EXPECT_NEAR( R(2,1) , 0.0 , 0.00001 ) ;
  EXPECT_NEAR( R(2,2) , 1.0 , 0.00001 ) ; 
}


TEST( icp , icp_rot_translate )
{
  // Generate point set 
  Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> target ;
  Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> data ; 

  // A sphere 
  const double pi = 3.14159265358979 ;
  const int nb_theta = 32 ;
  const int nb_phi   = 16 ;
  const double delta_theta = ( 2.0 * pi ) / nb_theta ;
  const double delta_phi = ( pi ) / nb_phi ;
  const double rad = 2.0 ;

  int nb_pt = 0 ;
  for( double theta = 0.0 ; theta < 2.0 * pi ; theta += delta_theta )
  {
    for( double phi = -pi/2.0 ; phi < pi/2.0 ; phi += delta_phi )
    {
      ++nb_pt ; 
    }
  }
  
  target.resize( nb_pt , 3 ) ;
  data.resize( nb_pt , 3 ) ; 

  openMVG::Mat3 r ;
  r = Eigen::AngleAxis<double>( pi / 3.0 , openMVG::Vec3( 1 , 1 , 1 ).normalized() ) ;

  const double tx = 10.0 ;
  const double ty = 12.0 ;
  const double tz = 3.0 ; 

  int id_pt = 0 ; 
  for( double theta = 0.0 ; theta < 2.0 * pi ; theta += delta_theta )
  {
    for( double phi = -pi/2.0 ; phi < pi/2.0 ; phi += delta_phi )
    {
      const double x = rad * std::cos( phi ) * std::cos( theta ) ;
      const double y = rad * std::cos( phi ) * std::sin( theta ) ;
      const double z = rad * std::sin( phi ) ;
    
      target( id_pt , 0 ) = x ;
      target( id_pt , 1 ) = y ;
      target( id_pt , 2 ) = z ; 


      const openMVG::Vec3 orig( x , y , z ) ;
      const openMVG::Vec3 tra = r * orig ; 
      
      data( id_pt , 0 ) = tra[0] + tx ; 
      data( id_pt , 1 ) = tra[1] + ty ;
      data( id_pt , 2 ) = tra[2] + tz ;  

      ++id_pt ;   
    }
  }

  // Inverse here 
  r = r.inverse().eval() ; 

  openMVG::Mat3 R ;
  openMVG::Vec3 t ; 

  openMVG::geometry::registration::ICP( target , data , 1000 , 0.0001 , t , R ) ; 

  std::cout << "Computed" << std::endl ; 
  std::cout << t << std::endl ; 
  std::cout << R << std::endl ; 

  std::cout << "Expected" << std::endl ; 
  std::cout << -tx << " - " << -ty << " - " << -tz << std::endl ; 
  std::cout << r << std::endl ; 

  EXPECT_NEAR( t[0] , -tx , 0.00001 ) ;
  EXPECT_NEAR( t[1] , -ty , 0.00001 ) ;
  EXPECT_NEAR( t[2] , -tz , 0.00001 ) ; 

  EXPECT_NEAR( R(0,0) , r(0,0) , 0.00001 ) ;
  EXPECT_NEAR( R(0,1) , r(0,1) , 0.00001 ) ;
  EXPECT_NEAR( R(0,2) , r(0,2) , 0.00001 ) ;

  EXPECT_NEAR( R(1,0) , r(1,0) , 0.00001 ) ;
  EXPECT_NEAR( R(1,1) , r(1,1) , 0.00001 ) ;
  EXPECT_NEAR( R(1,2) , r(1,2) , 0.00001 ) ;

  EXPECT_NEAR( R(2,0) , r(2,0) , 0.00001 ) ;
  EXPECT_NEAR( R(2,1) , r(2,1) , 0.00001 ) ;
  EXPECT_NEAR( R(2,2) , r(2,2) , 0.00001 ) ; 
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
