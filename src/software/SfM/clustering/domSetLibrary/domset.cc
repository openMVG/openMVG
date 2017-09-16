// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2016 nomoko AG, Srivathsan Murali<srivathsan@nomoko.camera>

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "domset.h"
#if OPENMVG_USE_OPENMP
#include <omp.h>
#endif

#include <fstream>
#include <limits>
#include <random>
#include <set>

#define for_parallel(i, nIters) for(int i = 0; i < static_cast<int>(nIters); ++i)

#include "openMVG/matching/matcher_kdtree_flann.hpp"
using namespace openMVG::matching;

namespace nomoko
{
void Domset::computeInformation()
{
  normalizePointCloud();
  voxelGridFilter( kVoxelSize, kVoxelSize, kVoxelSize );
  getAllDistances();
}

void Domset::normalizePointCloud()
{
  std::vector<Eigen::Vector3f> points_cpy;
  points_cpy.reserve(points.size());
  for (const auto & point_it : points)
    points_cpy.push_back(point_it.pos);

  ArrayMatcher_Kdtree_Flann<float> matcher;
  matcher.Build(points_cpy[0].data(), points_cpy.size(), 3);

  const size_t numResults( 1 );
  const size_t numPoints( points.size() );
  float totalDist = 0;
  pcCentre.pos << 0, 0, 0;
#if OPENMVG_USE_OPENMP
#pragma omp parallel for
#endif
  for_parallel( i, numPoints )
  {
    const Point &p           = points[ i ];
    const float queryPt[ 3 ] = {p.pos( 0 ), p.pos( 1 ), p.pos( 2 )};

    IndMatches ret_index( 2 );
    std::vector<float> out_dist_sqr( 2);
    const int NN = 2;
    matcher.SearchNeighbours(queryPt, 1, &ret_index, &out_dist_sqr, NN);

#if OPENMVG_USE_OPENMP
#pragma omp critical( TotalDistanceUpdate )
#endif
    {
      totalDist += ( p.pos - points[ ret_index[ 1 ].j_ ].pos ).norm();
      pcCentre.pos += ( p.pos / numPoints );
    }
  }

  // calculation the normalization scale
  const float avgDist = totalDist / numPoints;
  /*
  std::cerr << "Total distance = " << totalDist << std::endl;
  std::cerr << "Avg distance = " << avgDist << std::endl;
  */
  normScale = 1.f / avgDist;
  /*
  std::cerr << "Normalization Scale = " << normScale << std::endl;
  */

// normalizing the distances on points
#if OPENMVG_USE_OPENMP
#pragma omp parallel for
#endif
  for_parallel( i, numPoints )
  {
    points[ i ].pos = ( points[ i ].pos - pcCentre.pos ) * normScale;
  }

  // normalizing the camera center positions
  const size_t numViews( views.size() );
#if OPENMVG_USE_OPENMP
#pragma omp parallel for
#endif
  for_parallel( i, numViews )
  {
    views[ i ].trans = ( views[ i ].trans - pcCentre.pos ) * normScale;
  }
}

void Domset::deNormalizePointCloud()
{
  const size_t numPoints( points.size() );
// denormalizing points
#if OPENMVG_USE_OPENMP
#pragma omp parallel for
#endif
  for_parallel( i, numPoints )
  {
    points[ i ].pos = ( points[ i ].pos / normScale ) + pcCentre.pos;
  }

  const size_t numOldPoints( origPoints.size() );
#if OPENMVG_USE_OPENMP
#pragma omp parallel for
#endif
  for_parallel( i, numOldPoints )
  {
    origPoints[ i ].pos = ( origPoints[ i ].pos / normScale ) + pcCentre.pos;
  }

  // denormalizing camera centers
  const size_t numViews( views.size() );
#if OPENMVG_USE_OPENMP
#pragma omp parallel for
#endif
  for_parallel( i, numViews )
  {
    views[ i ].trans = ( views[ i ].trans / normScale ) + pcCentre.pos;
  }
}

void Domset::voxelGridFilter( const float &sizeX, const float &sizeY, const float &sizeZ )
{
  if ( sizeX <= 0.0f || sizeY <= 0.0f || sizeZ <= 0.0f )
  {
    std::cerr << "Invalid voxel grid dimensions error.\n";
    exit( 0 );
  }

  Point minPt;
  Point maxPt;
  const size_t numP = points.size();
  // finding the min and max values for the 3 dimensions
  const float mi = std::numeric_limits<float>::min();
  const float ma = std::numeric_limits<float>::max();
  minPt.pos << ma, ma, ma;
  maxPt.pos << mi, mi, mi;

  for ( size_t p = 0; p < numP; p++ )
  {
    const Point newSP = points[ p ];
    if ( newSP.pos( 0 ) < minPt.pos( 0 ) )
      minPt.pos( 0 ) = newSP.pos( 0 );
    if ( newSP.pos( 1 ) < minPt.pos( 1 ) )
      minPt.pos( 1 ) = newSP.pos( 1 );
    if ( newSP.pos( 2 ) < minPt.pos( 2 ) )
      minPt.pos( 2 ) = newSP.pos( 2 );
    if ( newSP.pos( 0 ) > maxPt.pos( 0 ) )
      maxPt.pos( 0 ) = newSP.pos( 0 );
    if ( newSP.pos( 1 ) > maxPt.pos( 1 ) )
      maxPt.pos( 1 ) = newSP.pos( 1 );
    if ( newSP.pos( 2 ) > maxPt.pos( 2 ) )
      maxPt.pos( 2 ) = newSP.pos( 2 );
  }

  // finding the number of voxels reqired
  //size_t numVoxelX = static_cast<size_t>( ceil( maxPt.pos( 0 ) - minPt.pos( 0 ) ) / sizeX );
  const size_t numVoxelY = static_cast<size_t>( ceil( maxPt.pos( 1 ) - minPt.pos( 1 ) ) / sizeY );
  const size_t numVoxelZ = static_cast<size_t>( ceil( maxPt.pos( 2 ) - minPt.pos( 2 ) ) / sizeZ );
  /*
  std::cout << "Max = " << maxPt.pos.transpose() << std::endl;
  std::cout << "Min = " << minPt.pos.transpose() << std::endl;
  std::cout << "Max - Min = " << ( maxPt.pos - minPt.pos ).transpose() << std::endl;
  std::cout << "VoxelSize X = " << sizeX << std::endl;
  std::cout << "VoxelSize Y = " << sizeY << std::endl;
  std::cout << "VoxelSize Z = " << sizeZ << std::endl;
  std::cout << "Number Voxel X = " << numVoxelX << std::endl;
  std::cout << "Number Voxel Y = " << numVoxelY << std::endl;
  std::cout << "Number Voxel Z = " << numVoxelZ << std::endl;
  */

  /* adding points to the voxels */
  std::map<size_t, std::vector<size_t>> voxels;
  std::vector<size_t> voxelIds;
#if OPENMVG_USE_OPENMP
#pragma omp parallel for
#endif
  for_parallel( p, numP )
  {
    const Point pt  = points[ p ];
    const size_t x  = static_cast<size_t>( floor( ( pt.pos( 0 ) - minPt.pos( 0 ) ) / sizeX ) );
    const size_t y  = static_cast<size_t>( floor( ( pt.pos( 1 ) - minPt.pos( 1 ) ) / sizeY ) );
    const size_t z  = static_cast<size_t>( floor( ( pt.pos( 2 ) - minPt.pos( 2 ) ) / sizeZ ) );
    const size_t id = ( z * numVoxelZ ) + ( y * numVoxelY ) + x;
#if OPENMVG_USE_OPENMP
#pragma omp critical( voxelGridUpdate )
#endif
    {
      if ( voxels.find( id ) == voxels.end() )
      {
        voxels[ id ] = std::vector<size_t>();
        voxelIds.push_back( id );
      }

      voxels[ id ].push_back( p );
    }
  }

  std::vector<Point> newPoints;
  const size_t numVoxelMaps = voxelIds.size();
#if OPENMVG_USE_OPENMP
#pragma omp parallel for
#endif
  for_parallel( vmId, numVoxelMaps )
  {
    const size_t vId  = voxelIds[ vmId ];
    const size_t nPts = voxels[ vId ].size();
    if ( nPts == 0 )
      continue;

    Eigen::Vector3f pos;
    std::set<size_t> vl;
    for ( const auto &p : voxels[ vId ] )
    {
      const Point pt = points[ p ];
      pos += pt.pos;
      const size_t numV = pt.viewList.size();
      for ( size_t j = 0; j < numV; j++ )
        vl.insert( pt.viewList[ j ] );
    }
    pos /= (float)nPts;

    Point newSP;
    newSP.pos      = pos;
    newSP.viewList = std::vector<size_t>( vl.begin(), vl.end() );
#if OPENMVG_USE_OPENMP
#pragma omp critical( pointsUpdate )
#endif
    {
      for ( const size_t viewID : vl )
      {
        views[ viewID ].viewPoints.push_back( newPoints.size() );
      }
      newPoints.push_back( newSP );
    }
  }

//  std::cerr << "New points = " << newPoints.size() << std::endl;
  origPoints.clear();
  points.swap( origPoints );
  points.swap( newPoints );
  // std::cerr << "Number of points = " << points.size() << std::endl;
} // voxelGridFilter

Eigen::MatrixXf Domset::getSimilarityMatrix( std::map<size_t, size_t> &xId2vId )
{
//  std::cout << "Generating Similarity Matrix " << std::endl;
  const size_t numC = xId2vId.size();
  const size_t numP = points.size();
  if ( numC == 0 || numP == 0 )
  {
    std::cerr << "Invalid Data\n";
    exit( 0 );
  }
  const float medianDist = getDistanceMedian( xId2vId );
//  std::cout << "Median dists = " << medianDist << std::endl;
  Eigen::MatrixXf simMat;
  simMat.resize( numC, numC );
#if OPENMVG_USE_OPENMP
#if _OPENMP > 200505 // collapse is only accessible from OpenMP 3.0
#pragma omp parallel for collapse( 2 )
#else
#pragma omp parallel for
#endif
#endif
  for_parallel( xId1, numC )
  {
    for_parallel( xId2, numC )
    {
      const size_t vId1 = xId2vId[ xId1 ];
      const size_t vId2 = xId2vId[ xId2 ];
      if ( vId1 == vId2 )
      {
        simMat( xId1, xId2 ) = 0;
      }
      else
      {
        const View v2   = views[ vId2 ];
        const View v1   = views[ vId1 ];
        const float sv  = computeViewSimilaity( v1, v2 );
        const float sd  = computeViewDistance( vId1, vId2, medianDist );
        const float sim = sv * sd;
        simMat( xId1, xId2 ) = sim;
      }
    }
  }
  return simMat;
} // getSimilarityMatrix

float Domset::computeViewDistance( const size_t &vId1, const size_t &vId2, const float &medianDist )
{
  if ( vId1 == vId2 )
    return 1.f;
  const float vd = viewDists( vId1, vId2 );
  const float dm = 1.f + exp( -( vd - medianDist ) / medianDist );
  return 1.f / dm;
}
float Domset::getDistanceMedian( const std::map<size_t, size_t> &xId2vId )
{
//  std::cout << "Finding Distance Median\n";

  if ( xId2vId.empty() )
  {
    std::cerr << "No Views initialized \n";
    exit( 0 );
  }

  const size_t numC = xId2vId.size();
  std::vector<float> dists;
  dists.reserve( numC * numC - numC );
  // float totalDist = 0;
  for ( size_t i = 0; i < numC; i++ )
  {
    const size_t v1 = xId2vId.at( i );
    for ( size_t j = 0; j < numC; j++ )
    {
      if ( i == j )
        continue;
      const size_t v2 = xId2vId.at( j );
      dists.push_back( viewDists( v1, v2 ) );
    }
  }
  std::sort( dists.begin(), dists.end() );
  return dists[ dists.size() / 2 ];
} // getDistanceMedian

void Domset::getAllDistances()
{
//  std::cout << "Finding View Distances\n";
  const size_t numC = views.size();
  if ( numC == 0 )
  {
    std::cerr << "No Views initialized \n";
    exit( 0 );
  }
  viewDists.resize( numC, numC );
  for ( size_t i = 0; i < numC; i++ )
  {
    const auto v1 = views[ i ];
    for ( size_t j = 0; j < numC; j++ )
    {
      const auto v2    = views[ j ];
      const float dist = ( v1.trans - v2.trans ).norm();
      viewDists( i, j ) = dist;
    }
  }
}

void Domset::findCommonPoints( const View &v1, const View &v2,
                               std::vector<size_t> &commonPoints )
{
  commonPoints.clear();
  const size_t numVP1 = v1.viewPoints.size();
  const size_t numVP2 = v2.viewPoints.size();
  const size_t minNum = std::min( numVP1, numVP2 );

  //std::sort(v1.viewPoints.begin(), v1.viewPoints.end());
  //std::sort(v2.viewPoints.begin(), v2.viewPoints.end());
  commonPoints.resize( minNum );

  const auto it = std::set_intersection( v1.viewPoints.begin(), v1.viewPoints.end(),
                                         v2.viewPoints.begin(), v2.viewPoints.end(), commonPoints.begin() );
  commonPoints.resize( it - commonPoints.begin() );
} // findCommonPoints

const float Domset::computeViewSimilaity( const View &v1, const View &v2 )
{
  std::vector<size_t> commonPoints;
  findCommonPoints( v1, v2, commonPoints );
  const size_t numCP = commonPoints.size();

  float w = 0.f;
#if OPENMVG_USE_OPENMP
#pragma omp parallel for
#endif
  for_parallel( p, numCP )
  {
    const auto pId = commonPoints[ p ];
    //for( const auto pId : commonPoints ){
    const Eigen::Vector3f c1 = (v1.trans - points[ pId ].pos).normalized();
    const Eigen::Vector3f c2 = (v2.trans - points[ pId ].pos).normalized();
    const float angle    = acos( c1.dot( c2 ) );
    const float expAngle = exp( -( angle * angle ) / kAngleSigma_2 );
//std::cerr << angle <<  " = " << expAngle << std::endl;
#if OPENMVG_USE_OPENMP
#pragma omp atomic
#endif
    w += expAngle;
  }
  const float ans = w / numCP;
  return ( ans != ans ) ? 0 : ans;
} // computeViewSimilaity

void Domset::computeClustersAP( std::map<size_t, size_t> &xId2vId,
                                std::vector<std::vector<size_t>> &clusters )
{
  const size_t numX = xId2vId.size();
  if ( numX == 0 )
  {
    std::cout << "Invalid map size\n";
    exit( 0 );
  }

  Eigen::MatrixXf S = getSimilarityMatrix( xId2vId );
  Eigen::MatrixXf R = Eigen::MatrixXf::Zero( numX, numX );
  Eigen::MatrixXf A = Eigen::MatrixXf::Zero( numX, numX );

  const float minFloat = std::numeric_limits<float>::lowest();
  for ( size_t m = 0; m < kNumIter; m++ )
  {
    // compute responsibilities
    Eigen::MatrixXf Rold = R;
    Eigen::MatrixXf AS   = A + S;
    std::vector<size_t> I( numX );
    Eigen::VectorXf Y( numX );
#if OPENMVG_USE_OPENMP
#pragma omp parallel for
#endif
    for ( size_t i = 0; i < numX; i++ )
    {
      Y( i ) = AS.row( i ).maxCoeff( &I[ i ] );
      AS( i, I[ i ] ) = minFloat;
    }

    std::vector<size_t> I2( numX );
    Eigen::VectorXf Y2( numX );
#if OPENMVG_USE_OPENMP
#pragma omp parallel for
#endif
    for ( size_t i = 0; i < numX; i++ )
    {
      Y2( i ) = AS.row( i ).maxCoeff( &I2[ i ] );
    }

    R = S - Y.replicate( 1, numX );
    for ( size_t i = 0; i < numX; i++ )
      R( i, I[ i ] ) = S( i, I[ i ] ) - Y2( i );
    R = ( ( 1 - lambda ) * R.array() ) + ( lambda * Rold.array() );

    // compute responsibilities
    Eigen::MatrixXf Aold = A;

    Eigen::MatrixXf Rp = ( R.array() > 0 ).select( R, 0 );
#if OPENMVG_USE_OPENMP
#pragma omp parallel for
#endif
    for ( size_t i = 0; i < numX; i++ )
    {
      Rp( i, i ) = R( i, i );
    }

    Eigen::VectorXf sumRp = Rp.colwise().sum();

    A                  = sumRp.transpose().replicate( numX, 1 ) - Rp;
    Eigen::VectorXf dA = A.diagonal();

    A = ( A.array() < 0 ).select( A, 0 );
#if OPENMVG_USE_OPENMP
#pragma omp parallel for
#endif
    for ( size_t i = 0; i < numX; i++ )
    {
      A( i, i ) = dA( i );
    }

    A = ( ( 1 - lambda ) * A.array() ) + ( lambda * Aold.array() );
  }

  Eigen::VectorXf E = ( A.diagonal() + R.diagonal() );
  E                 = ( E.array() > 0 ).select( E, 0 );

  // getting initial clusters
  std::set<size_t> centers;
  std::map<size_t, std::vector<size_t>> clMap;
  for ( size_t i = 0; i < numX; i++ )
  {
    if ( E( i ) > 0 )
    {
      centers.insert( i );
      clMap[ i ] = {};
    }
  }

  size_t idxForI = 0;
  for ( size_t i = 0; i < numX; i++ )
  {
    float maxSim = std::numeric_limits<float>::min();
    for ( auto const c : centers )
    {
      if ( S( i, c ) > maxSim )
      {
        idxForI = i;
        maxSim  = S( i, c );
      }
    }
    clMap[ idxForI ].push_back( i );
  }

  auto findCenter = [&]( const std::vector<size_t> &cl ) {
    std::vector<float> eValues;
    size_t i       = numX;
    float maxValve = -1;
    for ( size_t j : cl )
    {
      if ( maxValve < E( j ) )
      {
        maxValve = E( j );
        i        = j;
      }
    }
    return ( maxValve < 0 ) ? cl[ 0 ] : i;
  };
  auto addToNewClMap =
      []( std::map<size_t, std::vector<size_t>> &map, const std::vector<size_t> &cl, const size_t &center ) {
        if ( map.find( center ) == map.end() )
        {
          map[ center ] = cl;
        }
        else
        {
          map.at( center ).insert(
              map.at( center ).end(), cl.begin(), cl.end() );
        }
      };

  bool change = false;
  do
  {
    change = false;

    // enforcing min size constraints
    auto p1 = clMap.begin();
    while ( p1 != clMap.end() )
    {
      if ( p1->second.size() < kMinClusterSize )
      {
        float minDist     = std::numeric_limits<float>::max();
        int minId         = -1;
        const size_t vId1 = xId2vId.at( p1->first );
        for ( auto p2 = clMap.begin(); p2 != clMap.end(); ++p2 )
        {
          if ( p1->first == p2->first )
            continue;
          const size_t vId2 = xId2vId.at( p2->first );
          if ( viewDists( vId1, vId2 ) < minDist && ( p1->second.size() + p2->second.size() ) < kMaxClusterSize )
          {
            minDist = viewDists( vId1, vId2 );
            minId   = p2->first;
          }
        }
        if ( minId > -1 )
        {
          change = true;
          clMap[ minId ].insert( clMap[ minId ].end(),
                                 p1->second.begin(), p1->second.end() );
          p1 = clMap.erase( p1 );
        }
        else
        {
          ++p1;
        }
      }
      else
      {
        ++p1;
      }
    }

    // enforcing max size constraints
    for ( auto p = clMap.begin(); p != clMap.end(); ++p )
    {
      std::vector<size_t> cl = p->second;
      if ( cl.size() > kMaxClusterSize )
      {
        p         = clMap.erase( p );
        change    = true;
        auto it   = cl.begin();
        auto stop = it + kMaxClusterSize;
        while ( it < cl.end() )
        {
          auto tmp = std::vector<size_t>( it, stop );
          std::sort( tmp.begin(), tmp.end() );
          auto center = findCenter( tmp );
          addToNewClMap( clMap, tmp, center );
          it   = stop;
          stop = it + kMaxClusterSize;
          stop = ( stop > cl.end() ) ? cl.end() : stop;
        }
      }
    }
  } while ( change );

  // adding it to clusters vector
  for ( auto p = clMap.begin(); p != clMap.end(); ++p )
  {
    std::vector<size_t> cl;
    for ( const auto i : p->second )
    {
      cl.push_back( xId2vId[ i ] );
    }
    clusters.emplace_back( cl );
  }
}

void Domset::clusterViews( std::map<size_t, size_t> &xId2vId, const size_t &minClusterSize,
                           const size_t &maxClusterSize )
{
//  std::cout << "[ Clustering Views ] " << std::endl;
  const size_t umC = views.size();
  kMinClusterSize  = minClusterSize;
  kMaxClusterSize  = maxClusterSize;

  std::vector<std::vector<size_t>> clusters;
  computeClustersAP( xId2vId, clusters );

  deNormalizePointCloud();
  finalClusters.swap( clusters );
}

void Domset::clusterViews(
    const size_t &minClusterSize, const size_t &maxClusterSize )
{
//  std::cout << "[ Clustering Views ] " << std::endl;
  const size_t numC = views.size();
  kMinClusterSize   = minClusterSize;
  kMaxClusterSize   = maxClusterSize;

  std::map<size_t, size_t> xId2vId;
  for ( size_t i = 0; i < numC; i++ )
  {
    xId2vId[ i ] = i;
  }
  std::vector<std::vector<size_t>> clusters;
  computeClustersAP( xId2vId, clusters );

  deNormalizePointCloud();
  finalClusters.swap( clusters );
}

void Domset::printClusters()
{
  std::stringstream ss;
  ss << "Clusters : \n";
  for ( const auto cl : finalClusters )
  {
    ss << cl.size() << " : ";
    for ( const auto id : cl )
    {
      ss << id << " ";
    }
    ss << "\n\n";
  }
  std::cout << "Number of clusters = " << finalClusters.size() << std::endl;
  std::cout << ss.str();
}

void Domset::exportToPLY( const std::string &plyFilename, bool exportPoints )
{
  std::stringstream plys;
  plys << "ply\n"
       << "format ascii 1.0\n";

  size_t totalViews = 0;
  for ( const auto cl : finalClusters )
    totalViews += cl.size();
  const size_t numPts = origPoints.size();

  size_t totalPoints = totalViews;
  if ( exportPoints )
    totalPoints += numPts;
  plys << "element vertex "
       << totalPoints << std::endl
       << "property float x\n"
       << "property float y\n"
       << "property float z\n"
       << "property uchar red\n"
       << "property uchar green\n"
       << "property uchar blue\n"
       << "end_header\n";

  std::mt19937 gen;
  std::uniform_int_distribution<> dis(0, 255);
  for ( const auto cl : finalClusters )
  {
    const unsigned int
        red   = dis(gen),
        green = dis(gen),
        blue  = dis(gen);
    for ( const auto id : cl )
    {
      const auto & pos = views[ id ].trans;
      plys
          << pos( 0 ) << " " << pos( 1 ) << " " << pos( 2 ) << " "
          << red << " " << green << " " << blue << std::endl;
    }
  }

  if ( exportPoints )
  {
    for ( const auto pt : origPoints )
    {
      const auto &pos = pt.pos;

      plys << pos( 0 ) << " " << pos( 1 ) << " " << pos( 2 )
           << " 255 255 255" << std::endl;
    }
  }

  std::ofstream plyFile( plyFilename );
  if ( !plyFile.is_open() )
  {
    std::cout << "Cant open " << plyFilename << " file\n";
  }
  else
  {
    plyFile << plys.str();
    plyFile.close();
  }
}

}
