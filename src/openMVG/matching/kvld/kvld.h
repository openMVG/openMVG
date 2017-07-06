#ifndef OPENMVG_MATCHING_KVLD_H
#define OPENMVG_MATCHING_KVLD_H
/** @Main KVLD algorithm implementation
 ** @Containing scale image pyramid, VLD structure and KVLD algorithm
 ** @author Zhe Liu
 **/

/*
Copyright (C) 2011-12 Zhe Liu and Pierre Moulon.
All rights reserved.

This file is part of the KVLD library and is made available under
the terms of the BSD license (see the COPYING file).
*/

#include <algorithm>
#include <iostream>
#include <vector>

#include "openMVG/numeric/eigen_alias_definition.hpp"
#include "openMVG/matching/kvld/algorithm.h"

namespace openMVG { namespace image { template <typename T> class Image; } }

//Parameters concerning speed and performance
const bool uniqueMatch      = true;//if activated, a point can be matched to only one point in the other image. Note: if false, it also desactivate partially geometric verification
const double juge           = 0.35;
const size_t max_connection = 20;
const double distance_thres = 0.5;
const float min_dist        = 10;
const float maxContrast     = 300.0f;

//===inner parameters of VLD, usually not to change
const int dimension         = 10; //number of simplified SIFT-like chain used in a single vld
const int subdirection      = 8; // number of bins in a SIFT-like chain histogram
const int binNum            = 24;//number of bins for SIFT-like chain main direction. Must be a pair number

//===== initialize parameters for a KVLD process ====//
// inlierRate: the minimum rate down to which the KVLD should be reprocessed with a lower number of inlierRate, initially set to 0.04
// K: the minimum number of gvld-consistent (or vld-consistent) neighbors to select a match as correct one
// geometry: if true, KVLD will also take geometric verification into account. c.f. paper
//           if false, KVLD execute a pure photometric verification
struct KvldParameters
{
  float inlierRate;
  size_t K;
  bool geometry;
  KvldParameters(): inlierRate( 0.04 ), K( 3 ), geometry( true ){}
};

//====== Pyramid of scale images ======//
// elements as angles[0] and magnitudes[0] presents orientations and gradient norms of the original scale (biggest), and the following elements
// are scaled down by a step of sqrt(2.0) in size
//
// angles: store orientations of pixels of each scale image into a vector of images, which varies from 0 to 2*PI for each pixel
// magnitudes: store gradient norms of pixels of each scale image into a vector of images
struct ImageScale
{
  std::vector< openMVG::image::Image< float > > angles;
  std::vector< openMVG::image::Image< float > > magnitudes;
  std::vector< double > ratios;
  double radius_size;
  double step;

  ImageScale(const openMVG::image::Image< float >& I, double r = 5.0);
  int getIndex( const double r )const;

private:
  void GradAndNorm(
    const openMVG::image::Image< float >& I,
    openMVG::image::Image< float >& angle,
    openMVG::image::Image< float >& m);
};

//====== VLD structures ======//
class VLD
{
  double contrast;
  float distance;

  float begin_point[ 2 ];
  float end_point[ 2 ];

  Eigen::Matrix< int, dimension, 1 > principleAngle; //relative angle
  Eigen::Matrix< double, dimension, 1 > weight;
  Eigen::Matrix< double, dimension * subdirection, 1 > descriptor;//relative angle

public:
  inline double get_contrast()const{ return contrast; }
  //====================constructors=====================//
  template< typename T >
  VLD( const ImageScale& series, T const& P1, T const& P2 );
//=========================================class functions==============================================//
  inline double get_orientation()const
  {
    const float dy = end_point[ 1 ] - begin_point[ 1 ];
    const float  dx = end_point[ 0 ] - begin_point[ 0 ];
    float angle;
    if (anglefrom( dx, dy, angle ))
      return angle;
    else
      return 0.0;
  }
  inline double difference( const  VLD& vld2 )const
  {
    double diff[ 2 ];
    diff[ 0 ] = 0;
    diff[ 1 ] = 0;

    if (contrast > 300 || vld2.contrast > 300  || contrast <= 0 || vld2.contrast <=0 )
      return 128;

    for (int i = 0; i < dimension; i++ )
    {
      for (int j = 0; j < subdirection; j++ )
      {// term of descriptor
        diff[ 0 ] += std::abs( descriptor[ i * subdirection + j ] - vld2.descriptor[ i * subdirection + j ] );
      }
      //term of main SIFT like orientation
      diff[ 1 ] += std::min( std::abs( principleAngle[ i ] - vld2.principleAngle[ i ] ),
        binNum - std::abs( principleAngle[ i ] - vld2.principleAngle[ i ] ) ) * ( weight[ i ] + vld2.weight[ i ] );// orientation term
    }

    diff[ 0 ] *= 0.36;
    diff[ 1 ] *= 0.64 / ( binNum );
    //std::cout<<"diff = "<<diff[0]<<" "<<diff[1]<<std::endl;
    return diff[ 0 ] + diff[ 1 ];
  }

  inline void test() const
  {
    std::cout << "contrast = " << contrast << std::endl;
    std::cout << std::endl << "distance= " << distance << std::endl;

    std::cout << "weights   : ";
    for (int i = 0; i < dimension; i++ )
    {
      std::cout << weight[ i ] << " ";
    }
    std::cout << std::endl;
    for (int i = 0; i < dimension; i++ )
    {
      //cout<<"principle= "<<principleAngle[i]<<endl;
      for (int j = 0; j < subdirection; j++ )
      {
        std::cout << descriptor[ i * subdirection + j ] << " ";
      }
      std::cout << std::endl;
    }
  }

};

//==================KVLD algorithm======================//
//Output specification: out put matches are 1-to-1 matches, multiple matches with one point is not allowed (or you can deactive this verification part in the code)
//
//I1, I2: input images,
//
//F1, F2: list of keypoints, each keypoint should contain variables as: x, y (plus angle and scale if geometric verification is set to true in kvldParameters)
//
//matches: initial possible matches containing pairs of INDEX of matches features from F1 and F2. An element of matches as (1,3) present the 1st feature in F1
//           and the 3rd feature in F2
//
//matchesFiltered: the output list of matches filled by KVLD process.(it will be cleared to empty at the beginning of KVLD).
//
//score: list of score for output matches. If geometric verification is set to true, each element presents the geometric consistancy score of the corresponding match with its neighbor matches
//        otherwise, it presents the average photometric consistency score with its neighbor matches.
//
//E: gvld(or vld)-consistency matrix, for illustration reason, it has been externalized as an input of KVLD. it should be initialized to be a matches.size*matche.sizes table with all equal to -1
//   e.g.  Mat E = Mat::ones(matchesPair.size(),matchesPair.size())*(-1);
//
//valide: indices of whether the i th match in the initial match list is selected, for illustration reason, it has been externalized as an input of KVLD. it should be initialized to be a
//    matches.size vector with all equal to true.  e.g.  std::vector<bool> valide(size, true);
//
//kvldParameters: container of minimum inlier rate, the value of K (=3 initially) and geometric verification flag (true initially)

float KVLD(const openMVG::image::Image< float >& I1,
  const openMVG::image::Image< float >& I2,
  const std::vector<openMVG::features::SIOPointFeature> & F1,
  const std::vector<openMVG::features::SIOPointFeature> & F2,
  const std::vector< openMVG::Pair >& matches,
  std::vector< openMVG::Pair >& matchesFiltered,
  std::vector< double >& score,
  openMVG::Mat& E,
  std::vector< bool >& valide,
  KvldParameters& kvldParameters );

#endif // OPENMVG_MATCHING_KVLD_H
