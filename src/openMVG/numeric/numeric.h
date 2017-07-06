
// Copyright (c) 2007, 2008 libmv authors.
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_NUMERIC_NUMERIC_H
#define OPENMVG_NUMERIC_NUMERIC_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <numeric>
#include <string>
#include <vector>

#include "openMVG/numeric/eigen_alias_definition.hpp"

namespace openMVG
{

//--------------
//-- Function --
//--------------


/**
* @brief Compute square of a number
* @tparam T Type of the number to square
* @param x Input number
* @return square of x
*/
template<typename T>
inline T Square( T x )
{
  return x * x;
}


/**
* @brief Clamp value inside a given range
* @tparam T working type
* @param val Value to clamp
* @param min Lower bound of clamping range
* @param max Upper bound od clamping range
* @return clamped value
* @note Assuming range form a valid range (ie: min <= max )
*/
template<typename T>
inline T clamp( const T & val, const T& min, const T & max )
{
  return std::max( min, std::min( val, max ) );
  //(val < min) ? val : ((val>max) ? val : max);
}

/**
* @brief Given a vector, computes it's cross product matrix
*
* Cross product matrix is a helper matrix used to express cross product as a multiplication matrix \n
* Given two vectors \f$a=\begin{pmatrix}a_x\\a_y\\a_z\end{pmatrix}\f$ and \f$b=\begin{pmatrix}b_x\\b_y\\b_z\end{pmatrix}\f$, cross product \f$a\times b\f$ is equal to :\n
* \f$a\times b = \begin{pmatrix}a\end{pmatrix}_{\times} b = \begin{pmatrix}0 & -a_z & a_y \\ a_z & 0 & -a_x \\ -a_y & a_x & 0 \end{pmatrix} b \f$
*
* where \f$\begin{pmatrix}a\end{pmatrix}_{\times}\f$ is the cross product matrix of a.
*
* @param x Input vector
* @return Cross product matrix of a input vector
*/
Mat3 CrossProductMatrix( const Vec3 &x );


/**
* @brief Compute rotation matrix around X-axis
* @param angle Angle of rotation in radian
* @return Rotation matrix of given magnitude
*/
Mat3 RotationAroundX( double angle );


/**
* @brief Compute rotation matrix around Y-axis
* @param angle Angle of rotation in radian
* @return Rotation matrix of given magnitude
*/
Mat3 RotationAroundY( double angle );


/**
* @brief Compute rotation matrix around Z-axis
* @param angle Angle sof rotation in radian
* @return Rotation matrix of given magnitude
*/
Mat3 RotationAroundZ( double angle );


/**
* @brief Convert an angle from degree to radian
* @param degree Angle in degree
* @return Same angle in radian
* @note Assuming input angle is in range [0;360]
*/
inline double D2R( double degree )
{
  return degree * M_PI / 180.0;
}


/**
* @brief Convert an angle from radian to degree
* @param radian Angle in radian
* @return Same angle in degree
* @note Assuming input angle in range [0;2Pi]
*/
inline double R2D( double radian )
{
  return radian / M_PI * 180.0;
}

/**
* @brief Compute mean rotation magnitude of the given rotation matrix
* @param R2 Input rotation matrix
* @return magnitude of the rotation (in radian)
* @note Assuming R2 is a correct rotation matrix
* @note Mean is computed using the matrix column dot products to an Identity matrix
*/
double  getRotationMagnitude( const Mat3 & R2 );

/**
* @brief Gives an indication of the sign of the input value
* @param x Value to test
* @retval 1.0 if value is nul or positive
* @retval -1.0 if value is negative
*/
inline double SIGN( double x )
{
  return x < 0.0 ? -1.0 : 1.0;
}

/**
* @brief Compute L infinity norm
* \f$ \| v \|_{\infty} = \max ( |v_0| , |v_1| , \dots , |v_n| ) \f$
* @param x Input vector
* @return L infinity norm of input vector
*/
template<typename TVec>
inline double NormLInfinity( const TVec &x )
{
  return x.array().abs().maxCoeff();
}

/**
* @brief Compute L infinity distance between two vectors
* @param x first vector
* @param y second vector
* @return distance between input vectors using L infinity norm
*/
template<typename TVec>
inline double DistanceLInfinity( const TVec &x, const TVec &y )
{
  return NormLInfinity( x - y );
}

/**
* @brief Compute look at matrix
* Make a rotation matrix such that center becomes the direction of the
* positive z-axis, and y is oriented close to up by default.
* @param center New direction (z-axis)
* @param up Desired up vector (y-axis)
* @return Rotation matrix
*/
Mat3 LookAt( const Vec3 &center, const Vec3 & up = Vec3::UnitY() );


/**
* @brief Compute generic look at matrix
* @param eyePosition3D New center of rotation
* @param center3D Position where matrix look at (center3D-eyePosition3D forms the new z-axis)
* @param upVector3D Desired up vector (y-axis)
* @return Rotation matrix conforming the given parameters
*/
Mat3 LookAt2( const Vec3 &eyePosition3D,
              const Vec3 &center3D = Vec3::Zero(),
              const Vec3 &upVector3D = Vec3::UnitY() );

#define SUM_OR_DYNAMIC(x,y) (x==Eigen::Dynamic||y==Eigen::Dynamic)?Eigen::Dynamic:(x+y)

template<typename Derived1, typename Derived2>
struct hstack_return
{
  using Scalar =typename Derived1::Scalar;
  enum
  {
    RowsAtCompileTime = Derived1::RowsAtCompileTime,
    ColsAtCompileTime = SUM_OR_DYNAMIC( Derived1::ColsAtCompileTime, Derived2::ColsAtCompileTime ),
    Options = Derived1::Flags & Eigen::RowMajorBit ? Eigen::RowMajor : 0,
    MaxRowsAtCompileTime = Derived1::MaxRowsAtCompileTime,
    MaxColsAtCompileTime = SUM_OR_DYNAMIC( Derived1::MaxColsAtCompileTime, Derived2::MaxColsAtCompileTime )
  };
  using type =
    Eigen::Matrix<
      Scalar,
      RowsAtCompileTime,
      ColsAtCompileTime,
      Options,
      MaxRowsAtCompileTime,
      MaxColsAtCompileTime>;
};

template<typename Derived1, typename Derived2>
typename hstack_return<Derived1, Derived2>::type
HStack ( const Eigen::MatrixBase<Derived1>& lhs, const Eigen::MatrixBase<Derived2>& rhs )
{
  typename hstack_return<Derived1, Derived2>::type res;
  res.resize( lhs.rows(), lhs.cols() + rhs.cols() );
  res << lhs, rhs;
  return res;
}


template<typename Derived1, typename Derived2>
struct vstack_return
{
  using Scalar =  typename Derived1::Scalar;
  enum
  {
    RowsAtCompileTime = SUM_OR_DYNAMIC( Derived1::RowsAtCompileTime, Derived2::RowsAtCompileTime ),
    ColsAtCompileTime = Derived1::ColsAtCompileTime,
    Options = Derived1::Flags & Eigen::RowMajorBit ? Eigen::RowMajor : 0,
    MaxRowsAtCompileTime = SUM_OR_DYNAMIC( Derived1::MaxRowsAtCompileTime, Derived2::MaxRowsAtCompileTime ),
    MaxColsAtCompileTime = Derived1::MaxColsAtCompileTime
  };
  using type =
    Eigen::Matrix<
      Scalar,
      RowsAtCompileTime,
      ColsAtCompileTime,
      Options,
      MaxRowsAtCompileTime,
      MaxColsAtCompileTime>;
};

template<typename Derived1, typename Derived2>
typename vstack_return<Derived1, Derived2>::type
VStack ( const Eigen::MatrixBase<Derived1>& lhs, const Eigen::MatrixBase<Derived2>& rhs )
{
  typename vstack_return<Derived1, Derived2>::type res;
  res.resize( lhs.rows() + rhs.rows(), lhs.cols() );
  res << lhs, rhs;
  return res;
}
#undef SUM_OR_DYNAMIC

/**
* @brief Compute Frobenius norm
* \f$ \| A \|_2 = \sqrt{ \sum_{i=1}^n \sum_{j=0}^m a_{ij}^2 } \f$
* @param Input A input matrix
* @return Frobenius norm of given matrix
*/
template<typename TMat>
inline double FrobeniusNorm( const TMat &A )
{
  return A.norm();
}

/**
* @brief Compute distance between two matrices using Frobenius norm
* @param A first matrix
* @param B second matrix
* @return Distance between the input matrices given Frobenius norm
*/
template<typename TMat>
inline double FrobeniusDistance( const TMat &A, const TMat &B )
{
  return FrobeniusNorm( A - B );
}


/**
* @brief Compute similarity of matrices given cosine similarity mesure
* \f$ \cos( A , B ) = \frac{ A . B }{ \| A \|_2 \| B \|_2 } \f$
* @param a First matrix
* @param b Second matrix
* @return cosine similarity mesure between the input matrices
*/
template<class TMat>
double CosinusBetweenMatrices( const TMat &a, const TMat &b )
{
  return ( a.array() * b.array() ).sum() /
         FrobeniusNorm( a ) / FrobeniusNorm( b );
}

/**
* @brief Compute per row mean and variance
* @param A input matrix
* @param[out] mean_pointer a pointer to a vector where mean values are stored
* @param[out] variance_pointer a pointer to a vector where variance values are stored
* @note mean_pointer and variance_pointer vector may be resized to store all values
*/
void MeanAndVarianceAlongRows( const Mat &A,
                               Vec *mean_pointer,
                               Vec *variance_pointer );


/**
* @brief Export a matrix to a file in Text mode
* @param mat Matrix to export
* @param filename Path to the file where matrix will be written
* @param sPrefix Prefix before content of the matrix
* @retval true if export is correct
* @retval false if there was an error during export
*/
bool exportMatToTextFile( const Mat & mat, const std::string & filename,
                          const std::string & sPrefix = "A" );


/**
* @brief Test if value is a finite one
* @param val Input parameter
* @retval 0 if input is not a finite one
* @retval non-zero value if input is a finite one
*/
inline int is_finite( const double val )
{
#ifdef _MSC_VER
  return _finite( val );
#else
  return std::isfinite( val );
#endif
}


/**
* @brief Compute min, max, mean and median value a a given range
* @param begin start computation iterator
* @param end end computation iterator
* @param[out] min Minimum value of range
* @param[out] max Maximum value of range
* @param[out] mean Mean value of range
* @param[out] median Median value of range
* @return true if the statistical values can be estimated
*/
template <typename Type, typename DataInputIterator>
bool minMaxMeanMedian( DataInputIterator begin, DataInputIterator end,
                       Type & min, Type & max, Type & mean, Type & median )
{
  if (std::distance( begin, end ) < 1 )
  {
    return false;
  }

  std::vector<Type> vec_val( begin, end );
  std::sort( vec_val.begin(), vec_val.end() );
  min = vec_val[0];
  max = vec_val[vec_val.size() - 1];
  mean = std::accumulate( vec_val.begin(), vec_val.end(), Type( 0 ) )
    / static_cast<Type>( vec_val.size() );
  median = vec_val[vec_val.size() / 2];
  return true;
}


/**
* @brief Display to standard output min, max, mean and median value of input range
* @param begin start of range
* @param end end of range
*/
template <typename Type, typename DataInputIterator>
void minMaxMeanMedian( DataInputIterator begin, DataInputIterator end )
{
  Type min, max, mean, median;
  minMaxMeanMedian( begin, end, min, max, mean, median );
  std::cout << "\n"
            << "\t min: " << min << "\n"
            << "\t mean: " << mean << "\n"
            << "\t median: " << median << std::endl
            << "\t max: " << max << std::endl;
}

/**
 ** Split a range [ a; b [ into a set of n ranges :
 [ a; c1 [ U [ c1; c2 [ U ... U [ c(n-1); b [
  **
  Output range vector only store [ a , c1 , c2 , ... , b ]

 ** if input range can't be split (range [a;b[ size is less than nb_split, only return [a;b[ range
 **
 ** @param range_start Start of range to split
 ** @param range_end End of range to split
 ** @param nb_split Number of desired split
 ** @param d_range Output splitted range
 **/
template < typename T >
void SplitRange( const T range_start , const T range_end , const int nb_split ,
                 std::vector< T > & d_range )
{
  const T range_length = range_end - range_start;
  if (range_length < nb_split )
  {
    d_range.push_back( range_start );
    d_range.push_back( range_end );
  }
  else
  {
    const T delta_range = range_length / nb_split;

    d_range.push_back( range_start );
    for (int i = 1; i < nb_split; ++i )
    {
      d_range.push_back( range_start + i * delta_range );
    }
    d_range.push_back( range_end );
  }
}


} // namespace openMVG


#endif  // OPENMVG_NUMERIC_NUMERIC_H
