#include "openMVG/geometry/io/ply.hpp"

#include "testing/testing.h"

using Vec3   = openMVG::Vec3;
using Vec3uc = openMVG::Vec3uc;

TEST( ply, write_ascii )
{
  std::vector<Vec3> pts;
  const double start_x = -2.0;
  const double start_y = -3.0;
  const double start_z = -4.0;
  const double end_x   = 2.0;
  const double end_y   = 3.0;
  const double end_z   = 4.0;

  const int nb_split_x = 4;
  const int nb_split_y = 5;
  const int nb_split_z = 6;

  const double delta_x = ( end_x - start_x ) / static_cast<double>( nb_split_x );
  const double delta_y = ( end_y - start_y ) / static_cast<double>( nb_split_y );
  const double delta_z = ( end_z - start_z ) / static_cast<double>( nb_split_z );

  for ( int id_x = 0; id_x < nb_split_x; ++id_x )
  {
    for ( int id_y = 0; id_y < nb_split_y; ++id_y )
    {
      for ( int id_z = 0; id_z < nb_split_z; ++id_z )
      {
        const double x = start_x + id_x * delta_x;
        const double y = start_y + id_y * delta_y;
        const double z = start_z + id_z * delta_z;

        pts.push_back( Vec3( x, y, z ) );
      }
    }
  }

  bool res = openMVG::geometry::io::PLYWrite( pts, "cube_ascii.ply", false );
  EXPECT_EQ( res, true );
}

TEST( ply, write_ascii_col )
{
  std::vector<Vec3> pts;
  std::vector<Vec3uc> col;
  const double start_x = -2.0;
  const double start_y = -3.0;
  const double start_z = -4.0;
  const double end_x   = 2.0;
  const double end_y   = 3.0;
  const double end_z   = 4.0;

  const int nb_split_x = 4;
  const int nb_split_y = 5;
  const int nb_split_z = 6;

  const double delta_x = ( end_x - start_x ) / static_cast<double>( nb_split_x );
  const double delta_y = ( end_y - start_y ) / static_cast<double>( nb_split_y );
  const double delta_z = ( end_z - start_z ) / static_cast<double>( nb_split_z );

  for ( int id_x = 0; id_x < nb_split_x; ++id_x )
  {
    for ( int id_y = 0; id_y < nb_split_y; ++id_y )
    {
      for ( int id_z = 0; id_z < nb_split_z; ++id_z )
      {
        const double x = start_x + id_x * delta_x;
        const double y = start_y + id_y * delta_y;
        const double z = start_z + id_z * delta_z;

        pts.push_back( Vec3( x, y, z ) );
        col.push_back( Vec3uc( id_x * 255 / nb_split_x, id_y * 255 / nb_split_y, id_z * 255 / nb_split_z ) );
      }
    }
  }

  bool res = openMVG::geometry::io::PLYWrite( pts, nullptr, &col, "cube_ascii_col.ply", false );
  EXPECT_EQ( res, true );
}

TEST( ply, write_binary )
{
  std::vector<Vec3> pts;
  const double start_x = -2.0;
  const double start_y = -3.0;
  const double start_z = -4.0;
  const double end_x   = 2.0;
  const double end_y   = 3.0;
  const double end_z   = 4.0;

  const int nb_split_x = 4;
  const int nb_split_y = 5;
  const int nb_split_z = 6;

  const double delta_x = ( end_x - start_x ) / static_cast<double>( nb_split_x );
  const double delta_y = ( end_y - start_y ) / static_cast<double>( nb_split_y );
  const double delta_z = ( end_z - start_z ) / static_cast<double>( nb_split_z );

  for ( int id_x = 0; id_x < nb_split_x; ++id_x )
  {
    for ( int id_y = 0; id_y < nb_split_y; ++id_y )
    {
      for ( int id_z = 0; id_z < nb_split_z; ++id_z )
      {
        const double x = start_x + id_x * delta_x;
        const double y = start_y + id_y * delta_y;
        const double z = start_z + id_z * delta_z;

        pts.push_back( Vec3( x, y, z ) );
      }
    }
  }

  bool res = openMVG::geometry::io::PLYWrite( pts, "cube_bin.ply", true );
  EXPECT_EQ( res, true );
}

TEST( ply, write_binary_col )
{
  std::vector<Vec3> pts;
  std::vector<Vec3uc> col;
  const double start_x = -2.0;
  const double start_y = -3.0;
  const double start_z = -4.0;
  const double end_x   = 2.0;
  const double end_y   = 3.0;
  const double end_z   = 4.0;

  const int nb_split_x = 4;
  const int nb_split_y = 5;
  const int nb_split_z = 6;

  const double delta_x = ( end_x - start_x ) / static_cast<double>( nb_split_x );
  const double delta_y = ( end_y - start_y ) / static_cast<double>( nb_split_y );
  const double delta_z = ( end_z - start_z ) / static_cast<double>( nb_split_z );

  for ( int id_x = 0; id_x < nb_split_x; ++id_x )
  {
    for ( int id_y = 0; id_y < nb_split_y; ++id_y )
    {
      for ( int id_z = 0; id_z < nb_split_z; ++id_z )
      {
        const double x = start_x + id_x * delta_x;
        const double y = start_y + id_y * delta_y;
        const double z = start_z + id_z * delta_z;

        pts.push_back( Vec3( x, y, z ) );
        col.push_back( Vec3uc( id_x * 255 / nb_split_x, id_y * 255 / nb_split_y, id_z * 255 / nb_split_z ) );
      }
    }
  }

  bool res = openMVG::geometry::io::PLYWrite( pts, nullptr, &col, "cube_bin_col.ply", true );
  EXPECT_EQ( res, true );
}

TEST( ply, write_read_ascii )
{
  std::vector<Vec3> pts;
  const double start_x = -2.0;
  const double start_y = -3.0;
  const double start_z = -4.0;
  const double end_x   = 2.0;
  const double end_y   = 3.0;
  const double end_z   = 4.0;

  const int nb_split_x = 4;
  const int nb_split_y = 5;
  const int nb_split_z = 6;

  const double delta_x = ( end_x - start_x ) / static_cast<double>( nb_split_x );
  const double delta_y = ( end_y - start_y ) / static_cast<double>( nb_split_y );
  const double delta_z = ( end_z - start_z ) / static_cast<double>( nb_split_z );

  for ( int id_x = 0; id_x < nb_split_x; ++id_x )
  {
    for ( int id_y = 0; id_y < nb_split_y; ++id_y )
    {
      for ( int id_z = 0; id_z < nb_split_z; ++id_z )
      {
        const double x = start_x + id_x * delta_x;
        const double y = start_y + id_y * delta_y;
        const double z = start_z + id_z * delta_z;

        pts.push_back( Vec3( x, y, z ) );
      }
    }
  }

  bool res = openMVG::geometry::io::PLYWrite( pts, "cube_rw_ascii.ply", false );
  EXPECT_EQ( res, true );
  std::vector<Vec3> pts_r;
  res = openMVG::geometry::io::PLYRead( "cube_rw_ascii.ply", pts_r );
  EXPECT_EQ( res, true );
  EXPECT_EQ( pts.size(), pts_r.size() );

  for ( size_t id_point = 0; id_point < pts.size(); ++id_point )
  {
    EXPECT_NEAR( pts[ id_point ][0], pts_r[ id_point ][0] , 0.00001 );
    EXPECT_NEAR( pts[ id_point ][1], pts_r[ id_point ][1] , 0.00001 );
    EXPECT_NEAR( pts[ id_point ][2], pts_r[ id_point ][2] , 0.00001 );
  }
}

TEST( ply, write_read_ascii_col )
{
  std::vector<Vec3> pts;
  std::vector<Vec3uc> col;
  const double start_x = -2.0;
  const double start_y = -3.0;
  const double start_z = -4.0;
  const double end_x   = 2.0;
  const double end_y   = 3.0;
  const double end_z   = 4.0;

  const int nb_split_x = 4;
  const int nb_split_y = 5;
  const int nb_split_z = 6;

  const double delta_x = ( end_x - start_x ) / static_cast<double>( nb_split_x );
  const double delta_y = ( end_y - start_y ) / static_cast<double>( nb_split_y );
  const double delta_z = ( end_z - start_z ) / static_cast<double>( nb_split_z );

  for ( int id_x = 0; id_x < nb_split_x; ++id_x )
  {
    for ( int id_y = 0; id_y < nb_split_y; ++id_y )
    {
      for ( int id_z = 0; id_z < nb_split_z; ++id_z )
      {
        const double x = start_x + id_x * delta_x;
        const double y = start_y + id_y * delta_y;
        const double z = start_z + id_z * delta_z;

        pts.push_back( Vec3( x, y, z ) );
        col.push_back( Vec3uc( id_x * 255 / nb_split_x, id_y * 255 / nb_split_y, id_z * 255 / nb_split_z ) );
      }
    }
  }

  bool res = openMVG::geometry::io::PLYWrite( pts, nullptr, &col, "cube_rw_ascii_col.ply", false );
  EXPECT_EQ( res, true );
  std::vector<Vec3> pts_r;
  std::vector<Vec3uc> col_r; 
  res = openMVG::geometry::io::PLYRead( "cube_rw_ascii_col.ply", pts_r , nullptr , &col_r );
  EXPECT_EQ( res, true );
  EXPECT_EQ( pts.size(), pts_r.size() );
  EXPECT_EQ( col.size() , col_r.size() ) ; 

  for ( size_t id_point = 0; id_point < pts.size(); ++id_point )
  {
    EXPECT_NEAR( pts[ id_point ][0], pts_r[ id_point ][0] , 0.00001 );
    EXPECT_NEAR( pts[ id_point ][1], pts_r[ id_point ][1] , 0.00001 );
    EXPECT_NEAR( pts[ id_point ][2], pts_r[ id_point ][2] , 0.00001 );

    EXPECT_EQ( col[id_point] , col_r[id_point] ) ; 
  }
}

TEST( ply, write_read_ascii_col_nor )
{
  std::vector<Vec3> pts;
  std::vector<Vec3> nor;
  std::vector<Vec3uc> col;
  const double start_x = -2.0;
  const double start_y = -3.0;
  const double start_z = -4.0;
  const double end_x   = 2.0;
  const double end_y   = 3.0;
  const double end_z   = 4.0;

  const int nb_split_x = 4;
  const int nb_split_y = 5;
  const int nb_split_z = 6;

  const double delta_x = ( end_x - start_x ) / static_cast<double>( nb_split_x );
  const double delta_y = ( end_y - start_y ) / static_cast<double>( nb_split_y );
  const double delta_z = ( end_z - start_z ) / static_cast<double>( nb_split_z );

  for ( int id_x = 0; id_x < nb_split_x; ++id_x )
  {
    for ( int id_y = 0; id_y < nb_split_y; ++id_y )
    {
      for ( int id_z = 0; id_z < nb_split_z; ++id_z )
      {
        const double x = start_x + id_x * delta_x;
        const double y = start_y + id_y * delta_y;
        const double z = start_z + id_z * delta_z;

        pts.push_back( Vec3( x, y, z ) );
        col.push_back( Vec3uc( id_x * 255 / nb_split_x, id_y * 255 / nb_split_y, id_z * 255 / nb_split_z ) );
        nor.push_back( Vec3( x , y , z ).normalized() );
      }
    }
  }

  bool res = openMVG::geometry::io::PLYWrite( pts, &nor, &col, "cube_rw_ascii_col.ply", false );
  EXPECT_EQ( res, true );
  std::vector<Vec3> pts_r;
  std::vector<Vec3> nor_r; 
  std::vector<Vec3uc> col_r; 
  res = openMVG::geometry::io::PLYRead( "cube_rw_ascii_col.ply", pts_r , &nor_r , &col_r );
  EXPECT_EQ( res, true );
  EXPECT_EQ( pts.size(), pts_r.size() );
  EXPECT_EQ( col.size() , col_r.size() ) ; 

  for ( size_t id_point = 0; id_point < pts.size(); ++id_point )
  {
    EXPECT_NEAR( pts[ id_point ][0], pts_r[ id_point ][0] , 0.00001 );
    EXPECT_NEAR( pts[ id_point ][1], pts_r[ id_point ][1] , 0.00001 );
    EXPECT_NEAR( pts[ id_point ][2], pts_r[ id_point ][2] , 0.00001 );

    EXPECT_NEAR( nor[ id_point ][0] , nor_r[ id_point ][0] , 0.00001 ) ;
    EXPECT_NEAR( nor[ id_point ][1] , nor_r[ id_point ][1] , 0.00001 ) ;
    EXPECT_NEAR( nor[ id_point ][2] , nor_r[ id_point ][2] , 0.00001 ) ;
    

    EXPECT_EQ( col[id_point] , col_r[id_point] ) ; 
  }
}


TEST( ply, write_read_binary )
{
  std::vector<Vec3> pts;
  const double start_x = -2.0;
  const double start_y = -3.0;
  const double start_z = -4.0;
  const double end_x   = 2.0;
  const double end_y   = 3.0;
  const double end_z   = 4.0;

  const int nb_split_x = 4;
  const int nb_split_y = 5;
  const int nb_split_z = 6;

  const double delta_x = ( end_x - start_x ) / static_cast<double>( nb_split_x );
  const double delta_y = ( end_y - start_y ) / static_cast<double>( nb_split_y );
  const double delta_z = ( end_z - start_z ) / static_cast<double>( nb_split_z );

  for ( int id_x = 0; id_x < nb_split_x; ++id_x )
  {
    for ( int id_y = 0; id_y < nb_split_y; ++id_y )
    {
      for ( int id_z = 0; id_z < nb_split_z; ++id_z )
      {
        const double x = start_x + id_x * delta_x;
        const double y = start_y + id_y * delta_y;
        const double z = start_z + id_z * delta_z;

        pts.push_back( Vec3( x, y, z ) );
      }
    }
  }

  bool res = openMVG::geometry::io::PLYWrite( pts, "cube_wr_bin.ply" );
  EXPECT_EQ( res, true );
  std::vector<Vec3> pts_r;
  res = openMVG::geometry::io::PLYRead( "cube_wr_bin.ply", pts_r );
  EXPECT_EQ( res, true );
  EXPECT_EQ( pts.size(), pts_r.size() );

  for ( size_t id_point = 0; id_point < pts.size(); ++id_point )
  {
    EXPECT_EQ( pts[ id_point ], pts_r[ id_point ] );
  }
}

TEST( ply, write_read_binary_col )
{
  std::vector<Vec3> pts;
  std::vector<Vec3uc> col;
  const double start_x = -2.0;
  const double start_y = -3.0;
  const double start_z = -4.0;
  const double end_x   = 2.0;
  const double end_y   = 3.0;
  const double end_z   = 4.0;

  const int nb_split_x = 4;
  const int nb_split_y = 5;
  const int nb_split_z = 6;

  const double delta_x = ( end_x - start_x ) / static_cast<double>( nb_split_x );
  const double delta_y = ( end_y - start_y ) / static_cast<double>( nb_split_y );
  const double delta_z = ( end_z - start_z ) / static_cast<double>( nb_split_z );

  for ( int id_x = 0; id_x < nb_split_x; ++id_x )
  {
    for ( int id_y = 0; id_y < nb_split_y; ++id_y )
    {
      for ( int id_z = 0; id_z < nb_split_z; ++id_z )
      {
        const double x = start_x + id_x * delta_x;
        const double y = start_y + id_y * delta_y;
        const double z = start_z + id_z * delta_z;

        pts.push_back( Vec3( x, y, z ) );
        col.push_back( Vec3uc( id_x * 255 / nb_split_x, id_y * 255 / nb_split_y, id_z * 255 / nb_split_z ) );
      }
    }
  }

  bool res = openMVG::geometry::io::PLYWrite( pts, nullptr, &col, "cube_rw_bin_col.ply", true );
  EXPECT_EQ( res, true );
  std::vector<Vec3> pts_r;
  std::vector<Vec3uc> col_r;
  res = openMVG::geometry::io::PLYRead( "cube_rw_bin_col.ply", pts_r, nullptr, &col_r );
  EXPECT_EQ( res, true );
  EXPECT_EQ( pts.size(), pts_r.size() );
  EXPECT_EQ( col.size(), col_r.size() );
  EXPECT_EQ( pts.size(), col.size() );

  for ( size_t id_point = 0; id_point < pts.size(); ++id_point )
  {
    EXPECT_EQ( pts[ id_point ], pts_r[ id_point ] );
    EXPECT_EQ( col[ id_point ], col_r[ id_point ] );
  }
}

TEST( ply, write_read_binary_col_nor )
{
  std::vector<Vec3> pts;
  std::vector<Vec3> nor; 
  std::vector<Vec3uc> col;
  const double start_x = -2.0;
  const double start_y = -3.0;
  const double start_z = -4.0;
  const double end_x   = 2.0;
  const double end_y   = 3.0;
  const double end_z   = 4.0;

  const int nb_split_x = 4;
  const int nb_split_y = 5;
  const int nb_split_z = 6;

  const double delta_x = ( end_x - start_x ) / static_cast<double>( nb_split_x );
  const double delta_y = ( end_y - start_y ) / static_cast<double>( nb_split_y );
  const double delta_z = ( end_z - start_z ) / static_cast<double>( nb_split_z );

  for ( int id_x = 0; id_x < nb_split_x; ++id_x )
  {
    for ( int id_y = 0; id_y < nb_split_y; ++id_y )
    {
      for ( int id_z = 0; id_z < nb_split_z; ++id_z )
      {
        const double x = start_x + id_x * delta_x;
        const double y = start_y + id_y * delta_y;
        const double z = start_z + id_z * delta_z;

        pts.push_back( Vec3( x, y, z ) );
        col.push_back( Vec3uc( id_x * 255 / nb_split_x, id_y * 255 / nb_split_y, id_z * 255 / nb_split_z ) );
        nor.push_back( Vec3(x,y,z).normalized() ) ;
      }
    }
  }

  bool res = openMVG::geometry::io::PLYWrite( pts, &nor , &col, "cube_rw_bin_col.ply", true );
  EXPECT_EQ( res, true );
  std::vector<Vec3> pts_r;
  std::vector<Vec3> nor_r; 
  std::vector<Vec3uc> col_r;
  res = openMVG::geometry::io::PLYRead( "cube_rw_bin_col.ply", pts_r, &nor_r , &col_r );
  EXPECT_EQ( res, true );
  EXPECT_EQ( pts.size(), pts_r.size() );
  EXPECT_EQ( col.size(), col_r.size() );
  EXPECT_EQ( pts.size(), col.size() );

  for ( size_t id_point = 0; id_point < pts.size(); ++id_point )
  {
    EXPECT_EQ( pts[ id_point ], pts_r[ id_point ] );
    EXPECT_EQ( nor[ id_point ], nor_r[ id_point ] );
    EXPECT_EQ( col[ id_point ], col_r[ id_point ] );
  }
}


/* ************************************************************************* */
int main()
{
  TestResult tr;
  return TestRegistry::runAllTests( tr );
}
/* ************************************************************************* */
