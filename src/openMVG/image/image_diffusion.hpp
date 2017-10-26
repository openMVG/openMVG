// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2014 Pierre Moulon, Romuald Perrot.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_IMAGE_IMAGE_DIFFUSION_HPP
#define OPENMVG_IMAGE_IMAGE_DIFFUSION_HPP

#ifdef _MSC_VER
  #pragma warning(once:4244)
#endif

#include <algorithm>
#include <vector>
#include <cmath>

#ifdef OPENMVG_USE_OPENMP
#include <omp.h>
#endif

#include "openMVG/numeric/numeric.h"

namespace openMVG
{
namespace image
{

/**
 ** Compute Perona and Malik G2 diffusion coefficient
 ** @param Lx Image of X-derivative
 ** @param Ly Image of Y-derivative
 ** @param k sensitivity factor
 ** @param out output coefficient
 ** NOTE : assume Lx and Ly have same size and image are in float format
 **/
template <typename Image>
void ImagePeronaMalikG2DiffusionCoef( const Image & Lx , const Image & Ly , const typename Image::Tpixel k , Image & out )
{
  const int width = Lx.Width();
  const int height = Lx.Height();

  if (width != out.Width() || height != out.Height())
  {
    out.resize( width , height );
  }

  using Real = typename Image::Tpixel;
  out.array() = ( static_cast<Real>( 1.f ) + ( Lx.array().square() + Ly.array().square() ) / ( k * k ) ).inverse();
}

/**
** Apply Fast Explicit Diffusion to an Image (on central part)
** @param src input image
** @param diff diffusion coefficient image
** @param half_t Half diffusion time
** @param out Output image
** @param row_start Row range beginning (range is [row_start; row_end [ )
** @param row_end Row range end (range is [row_start; row_end [ )
**/
template<typename Image>
void ImageFEDCentral( const Image & src , const Image & diff , const typename Image::Tpixel half_t , Image & out ,
                      const int row_start , const int row_end )
{
  using Real = typename Image::Tpixel;
  const int width = src.Width();
  Real n_diff[4];
  Real n_src[4];
  // Compute FED step on general range
  for (int i = row_start; i < row_end; ++i)
  {
    for (int j = 1; j < width - 1; ++j)
    {
      // Retrieve neighbors : TODO check if we need a cache efficient version ?
      n_diff[0] = diff( i , j + 1 );
      n_diff[1] = diff( i - 1 , j );
      n_diff[2] = diff( i , j - 1 );
      n_diff[3] = diff( i + 1 , j );
      n_src[0] = src( i , j + 1 );
      n_src[1] = src( i - 1 , j );
      n_src[2] = src( i , j - 1 );
      n_src[3] = src( i + 1 , j );

      // Compute diffusion factor for given pixel
      const Real cur_src = src( i , j );
      const Real cur_diff = diff( i , j );
      const Real a = ( cur_diff + n_diff[0] ) * ( n_src[0] - cur_src );
      const Real b = ( cur_diff + n_diff[1] ) * ( cur_src - n_src[1] );
      const Real c = ( cur_diff + n_diff[2] ) * ( cur_src - n_src[2] );
      const Real d = ( cur_diff + n_diff[3] ) * ( n_src[3] - cur_src );
      const Real value = half_t * ( a - c + d - b );
      out( i , j ) = value;
    }
  }
}

/**
** Apply Fast Explicit Diffusion to an Image (on central part)
** @param src input image
** @param diff diffusion coefficient image
** @param half_t Half diffusion time
** @param out Output image
**/
template<typename Image>
void ImageFEDCentralCPPThread( const Image & src , const Image & diff , const typename Image::Tpixel half_t , Image & out )
{
#ifdef OPENMVG_USE_OPENMP
  const int nb_thread = omp_get_max_threads();
#else
  const int nb_thread = 1;
#endif

  // Compute ranges
  std::vector<int > range;
  SplitRange( 1 , ( int ) ( src.rows() - 1 ) , nb_thread , range );

#ifdef OPENMVG_USE_OPENMP
  #pragma omp parallel for schedule(dynamic)
#endif
  for (int i = 1; i < static_cast<int>( range.size() ); ++i)
  {
    ImageFEDCentral( src, diff, half_t, out, range[i - 1] , range[i] );
  }
}

/**
** Apply Fast Explicit Diffusion of an Image
** @param src input image
** @param diff diffusion coefficient image
** @param t diffusion time
** @param out output image
**/
template<typename Image>
void ImageFED( const Image & src , const Image & diff , const typename Image::Tpixel t , Image & out )
{
  using Real = typename Image::Tpixel;
  const int width = src.Width();
  const int height = src.Height();
  const Real half_t = t * static_cast<Real>( 0.5 );
  if (out.Width() != width || out.Height() != height)
  {
    out.resize( width , height );
  }
  Real n_diff[4];
  Real n_src[4];

  // Take care of the central part
  ImageFEDCentralCPPThread( src , diff , half_t , out );

  // Take care of the border
  // - first/last row
  // - first/last col

  // Compute FED step on first row
  for (int j = 1; j < width - 1; ++j)
  {
    n_diff[0] = diff( 0 , j + 1 );
    n_diff[2] = diff( 0 , j - 1 );
    n_diff[3] = diff( 1 , j );
    n_src[0] = src( 0 , j + 1 );
    n_src[2] = src( 0 , j - 1 );
    n_src[3] = src( 1 , j );

    // Compute diffusion factor for given pixel
    const Real cur_src = src( 0 , j );
    const Real cur_diff = diff( 0 , j );
    const Real a = ( cur_diff + n_diff[0] ) * ( n_src[0] - cur_src );
    const Real c = ( cur_diff + n_diff[2] ) * ( cur_src - n_src[2] );
    const Real d = ( cur_diff + n_diff[3] ) * ( n_src[3] - cur_src );
    const Real value = half_t * ( a - c + d );
    out( 0 , j ) = value;
  }

  // Compute FED step on last row
  for (int j = 1; j < width - 1; ++j)
  {
    n_diff[0] = diff( height - 1 , j + 1 );
    n_diff[1] = diff( height - 2 , j );
    n_diff[2] = diff( height - 1 , j - 1 );
    n_src[0] = src( height - 1 , j + 1 );
    n_src[1] = src( height - 2 , j );
    n_src[2] = src( height - 1 , j - 1 );

    // Compute diffusion factor for given pixel
    const Real cur_src = src( height - 1 , j );
    const Real cur_diff = diff( height - 1 , j );
    const Real a = ( cur_diff + n_diff[0] ) * ( n_src[0] - cur_src );
    const Real b = ( cur_diff + n_diff[1] ) * ( cur_src - n_src[1] );
    const Real c = ( cur_diff + n_diff[2] ) * ( cur_src - n_src[2] );
    const Real value = half_t * ( a - c - b );
    out( height - 1 , j ) = value;
  }

  // Compute FED step on first col
  for (int i = 1; i < height - 1; ++i)
  {
    n_diff[0] = diff( i , 1 );
    n_diff[1] = diff( i - 1 , 0 );
    n_diff[3] = diff( i + 1 , 0 );
    n_src[0] = src( i , 1 );
    n_src[1] = src( i - 1 , 0 );
    n_src[3] = src( i + 1 , 0 );

    // Compute diffusion factor for given pixel
    const Real cur_src = src( i , 0 );
    const Real cur_diff = diff( i , 0 );
    const Real a = ( cur_diff + n_diff[0] ) * ( n_src[0] - cur_src );
    const Real b = ( cur_diff + n_diff[1] ) * ( cur_src - n_src[1] );
    const Real d = ( cur_diff + n_diff[3] ) * ( n_src[3] - cur_src );
    const Real value = half_t * ( a + d - b );
    out( i , 0 ) = value;
  }

  // Compute FED step on last col
  for (int i = 1; i < height - 1; ++i)
  {
    n_diff[1] = diff( i - 1 , width - 1 );
    n_diff[2] = diff( i , width - 2 );
    n_diff[3] = diff( i + 1 , width - 1 );
    n_src[1] = src( i - 1 , width - 1 );
    n_src[2] = src( i , width - 2 );
    n_src[3] = src( i + 1 , width - 1 );

    // Compute diffusion factor for given pixel
    const Real cur_src = src( i , width - 1 );
    const Real cur_diff = diff( i , width - 1 );
    const Real b = ( cur_diff + n_diff[1] ) * ( cur_src - n_src[1] );
    const Real c = ( cur_diff + n_diff[2] ) * ( cur_src - n_src[2] );
    const Real d = ( cur_diff + n_diff[3] ) * ( n_src[3] - cur_src );
    const Real value = half_t * ( - c + d - b );
    out( i , width - 1 ) = value;
  }
}

/**
 ** Compute Fast Explicit Diffusion cycle
 ** @param self input/output image
 ** @param diff diffusion coefficient
 ** @param tau cycle timing vector
 **/
template<typename Image>
void ImageFEDCycle( Image & self , const Image & diff , const std::vector<typename Image::Tpixel > & tau )
{
  Image tmp;
  for (int i = 0; i < tau.size(); ++i)
  {
    ImageFED( self , diff , tau[i] , tmp );
    self.array() += tmp.array();
  }
}

/**
* Compute if a number is prime of not
* @param i Input number to test
* @retval true if number is prime
* @retval false if number is not prime
* @todo Move this function elsewhere since it's not an image related function
*/
inline bool IsPrime( const int i )
{
  if (i == 1)
  {
    return false;
  }
  if (i == 2 || i == 3)
  {
    return true;
  }
  if (i % 2 == 0)
  {
    return false;
  }

  const size_t i_root = static_cast<int>( sqrt( static_cast<double>( i + 1 ) ) );

  for (size_t cur = 3; cur <= i_root; cur += 2)
  {
    if (i % cur == 0)
    {
      return false;
    }
  }
  return true;
}

/**
* @brief Get the next prime number greater or equal to input
* @param i Input number
* @return next prime greater or equal to input
*/
inline int NextPrimeGreaterOrEqualTo( const int i )
{
  if (IsPrime( i ))
  {
    return i;
  }
  else
  {
    int cur = i + 1;

    while (!IsPrime( cur ))
    {
      ++cur;
    }
    return cur;
  }
}

/**
 ** Compute FED cycle timings using total time
 ** @param T total time
 ** @param Tmax cycle stability limit (max : 0.25)
 ** @param tau vector of FED cycle timings
 ** @return number of cycle timings
 **/
template<typename Real>
int FEDCycleTimings( const Real T , const Real Tmax , std::vector<Real > & tau )
{
  // Number of timings
  const int n = ceil( sqrt( ( 3.0 * static_cast<double>( T ) ) / Tmax +  0.25 ) - 0.5 ) + 0.5;

  // Scaling factor
  const Real scale = 3.0 * T / ( Tmax * static_cast<Real>( n * n + n ) );

  // only constants
  const Real cos_fact = 1.0 / ( static_cast<Real>( 4 * n ) + 2.0 );
  const Real glo_fact = scale * Tmax / 2.0;

  // Compute cycle timings
  tau.resize( n );
  for (int j = 0; j < n; ++j)
  {
    const Real cos_j = cos( M_PI * ( static_cast<Real>( 2 * j + 1 ) ) * cos_fact );
    tau[ j ] = glo_fact / ( cos_j * cos_j );
  }

  // Compute Kappa reordering using kappa = n / 2
  const int kappa = n / 2;

  const int p = NextPrimeGreaterOrEqualTo( n + 1 );

  // Store new positions
  std::vector<Real > tmp( n );
  for (int i = 0 , k = 0; i < n; ++i , ++k)
  {
    // Search new index
    int index = n;
    while (( index = ( ( k + 1 ) * kappa ) % p - 1 ) >= n)
    {
      ++k;
    }

    tmp[ i ] = tau[ index ];
  }

  // Get new vector
  std::swap( tmp , tau );
  return n;
}

}  // namespace image
}  // namespace openMVG

#endif //  OPENMVG_IMAGE_IMAGE_DIFFUSION_HPP
