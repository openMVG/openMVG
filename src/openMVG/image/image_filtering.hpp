// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2014 Romuald Perrot.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_IMAGE_IMAGE_FILTERING_HPP
#define OPENMVG_IMAGE_IMAGE_FILTERING_HPP

/**
 ** @file
 ** standard image filtering functions :
 ** - X,Y derivatives using : Central difference, Sobel filter, Scharr filter
 ** - Gaussian blur
 **/

//------------------
//-- Bibliography --
//------------------
//- [1] "A Scheme for Coherence-Enhancing Diffusion Filtering with Optimized Rotation Invariance."
//- Authors: Joachim Weickert and Hanno Scharr.
//- Date: September 2002.
//- Journal : Journal of Visual Communication and Image Representation.

#include "openMVG/image/image_convolution.hpp"

namespace openMVG
{
namespace image
{

/**
 ** Compute X-derivative using central difference
 ** @param img Input image
 ** @param out Output image
 ** @param normalize true if kernel must be scaled by 1/2
 **/
template< typename Image >
void ImageXDerivative( const Image & img , Image & out , const bool normalize = true )
{
  Vec3 kernel( -1.0, 0.0, 1.0 );

  if (normalize )
  {
    kernel *= 0.5;
  }

  ImageHorizontalConvolution( img , kernel , out );
}


/**
 ** Compute Y-derivative using central difference
 ** @param img Input image
 ** @param out Output image
 ** @param normalize true if kernel must be normalized
 **/
template< typename Image >
void ImageYDerivative( const Image & img , Image & out , const bool normalize = true )
{
  Vec3 kernel( -1.0, 0.0, 1.0 );

  if (normalize )
  {
    kernel *= 0.5;
  }

  ImageVerticalConvolution( img , kernel , out );
}


/**
 ** Compute X-derivative using 3x3 Sobel kernel
 ** @param img Input image
 ** @param out Output image
 ** @param normalize true if kernel must be scaled by 1/8
 **/
template< typename Image >
void ImageSobelXDerivative( const Image & img , Image & out , const bool normalize = true )
{
  Vec3 kernel_horiz( -1.0, 0.0, 1.0 );

  if (normalize )
  {
    kernel_horiz *= 0.5;
  }

  Vec3 kernel_vert( 1.0, 2.0, 1.0 );

  if (normalize )
  {
    kernel_vert *= 0.25;
  }

  ImageSeparableConvolution( img , kernel_horiz , kernel_vert , out );
}


/**
 ** Compute Y-derivative using 3x3 Sobel kernel
 ** @param img Input image
 ** @param out Output image
 ** @param normalize true if kernel must be scaled by 1/8
 **/
template< typename Image >
void ImageSobelYDerivative( const Image & img , Image & out , const bool normalize = true )
{
  Vec3 kernel_horiz( 1.0, 2.0, 1.0 );

  if (normalize )
  {
    kernel_horiz *= 0.25;
  }

  Vec3 kernel_vert( -1.0, 0.0, 1.0 );

  if (normalize )
  {
    kernel_vert *= 0.5;
  }

  ImageSeparableConvolution( img , kernel_horiz , kernel_vert , out );
}

/**
 ** Compute X-derivative using 3x3 Scharr kernel
 ** @param img Input image
 ** @param out Output image
 ** @param normalize true if kernel must be scaled by 1/32
 **/
template< typename Image >
void ImageScharrXDerivative( const Image & img , Image & out , const bool normalize = true )
{
  Vec3 kernel_horiz( -1.0, 0.0, 1.0 );

  if (normalize )
  {
    kernel_horiz *= 0.5;
  }

  Vec3 kernel_vert( 3.0, 10.0, 3.0 );

  if (normalize )
  {
    kernel_vert *= 1.0 / 16.0;
  }

  ImageSeparableConvolution( img , kernel_horiz , kernel_vert , out );
}

/**
 ** Compute Y-derivative using 3x3 Scharr filter
 ** @param img Input image
 ** @param out Output image
 ** @param normalize true if kernel must be scaled by 1/32
 **/
template< typename Image >
void ImageScharrYDerivative( const Image & img , Image & out , const bool normalize = true )
{
  Vec3 kernel_horiz( 3.0, 10.0, 3.0 );

  if (normalize )
  {
    kernel_horiz *= 1.0 / 16.0;
  }

  Vec3 kernel_vert( -1.0, 0.0, 1.0 );

  if (normalize )
  {
    kernel_vert *= 0.5;
  }

  ImageSeparableConvolution( img , kernel_horiz , kernel_vert , out );
}


/**
 ** Compute X-derivative using scaled Scharr filter
 ** @param img Input image
 ** @param out Output image
 ** @param scale scale of filter (1 -> 3x3 filter; 2 -> 5x5, ...)
 ** @param bNormalize true if kernel must be normalized
 **/
template< typename Image >
void ImageScaledScharrXDerivative( const Image & img , Image & out , const int scale , const bool bNormalize = true )
{
  const int kernel_size = 3 + 2 * ( scale - 1 );

  Vec kernel_vert( kernel_size );
  Vec kernel_horiz( kernel_size );

  /*
  General X-derivative function
                              | -1   0   1 |
  D = 1 / ( 2 h * ( w + 2 ) ) | -w   0   w |
                              | -1   0   1 |
  */

  kernel_horiz.fill( 0.0 );
  kernel_horiz( 0 )               = -1.0;
  // kernel_horiz( kernel_size / 2 ) = 0.0;
  kernel_horiz( kernel_size - 1 ) = 1.0;

  // Scharr parameter for derivative
  const double w = 10.0 / 3.0;

  kernel_vert.fill( 0.0 );
  kernel_vert( 0 )               = 1.0;
  kernel_vert( kernel_size / 2 ) = w;
  kernel_vert( kernel_size - 1 ) = 1.0;

  if (bNormalize )
  {
    kernel_vert *= 1.0 / ( 2.0 * scale * ( w + 2.0 ) );
  }

  ImageSeparableConvolution( img , kernel_horiz , kernel_vert , out );
}



/**
 ** Compute Y-derivative using scaled Scharr filter
 ** @param img Input image
 ** @param out Output image
 ** @param scale scale of filter (1 -> 3x3 filter; 2 -> 5x5, ...)
 ** @param bNormalize true if kernel must be normalized
 **/
template< typename Image >
void ImageScaledScharrYDerivative( const Image & img , Image & out , const int scale , const bool bNormalize = true )
{
  /*
  General Y-derivative function
                              | -1  -w  -1 |
  D = 1 / ( 2 h * ( w + 2 ) ) |  0   0   0 |
                              |  1   w   1 |

  */
  const int kernel_size = 3 + 2 * ( scale - 1 );

  Vec kernel_vert( kernel_size );
  Vec kernel_horiz( kernel_size );

  // Scharr parameter for derivative
  const double w = 10.0 / 3.0;


  kernel_horiz.fill( 0.0 );
  kernel_horiz( 0 )               = 1.0;
  kernel_horiz( kernel_size / 2 ) = w;
  kernel_horiz( kernel_size - 1 ) = 1.0;

  if (bNormalize )
  {
    kernel_horiz *= 1.0 / ( 2.0 * scale * ( w + 2.0 ) );
  }

  kernel_vert.fill( 0.0 );
  kernel_vert( 0 ) = - 1.0;
  // kernel_vert( kernel_size / 2 ) = 0.0;
  kernel_vert( kernel_size - 1 ) = 1.0;

  ImageSeparableConvolution( img , kernel_horiz , kernel_vert , out );
}


/**
 ** Compute (isotropic) gaussian filtering of an image using filter width of k * sigma
 ** @param img Input image
 ** @param sigma standard deviation of kernel
 ** @param out Output image
 ** @param k confidence interval param - kernel is width k * sigma * 2 + 1 -- using k = 3 gives 99% of gaussian curve
 ** @param border_mgmt either BORDER_COPY or BORDER_CROP to tell what to do with borders
 **/
template< typename Image >
void ImageGaussianFilter( const Image & img , const double sigma , Image & out , const int k = 3 )
{
  // Compute Gaussian filter
  const int k_size    = ( int ) 2 * k * sigma + 1;
  const int half_k_size = k_size / 2;

  const double exp_scale = 1.0 / ( 2.0 * sigma * sigma );

  // Compute 1D Gaussian filter
  Vec kernel_horiz( k_size );

  double sum = 0;
  for (int i = 0; i < k_size; ++i )
  {
    const double dx = ( i - half_k_size );
    kernel_horiz( i ) = exp( - dx * dx * exp_scale );
    sum += kernel_horiz( i );
  }

  // Normalize kernel (to have \sum_i kernel_horiz( i ) = 1 and avoid energy loss)
  const double inv = 1.0 / sum;
  for (int i = 0; i < k_size; ++i )
  {
    kernel_horiz( i ) *= inv;
  }

  // Vertical kernel is the same as the horizontal one
  const Vec & kernel_vert = kernel_horiz;

  ImageSeparableConvolution( img , kernel_horiz , kernel_vert , out );
}

/**
 ** @brief Compute 1D gaussian kernel of specified width
 ** @param size Size of kernel (0 for automatic window)
 ** @param sigma Gaussian scale
 ** @return Kernel using specified parameters
 **/
inline Vec ComputeGaussianKernel( const size_t size , const double sigma )
{
  // If kernel size is 0 computes it's size using uber formula
  size_t k_size = ( size == 0 ) ? ceil( 2.0 * ( 1.0 + ( sigma - 0.8 ) / ( 0.3 ) ) ) : size;

  // Make kernel odd width
  k_size = ( k_size % 2 == 0 ) ? k_size + 1 : k_size;
  const size_t half_k_size = ( k_size - 1 ) / 2;

  Vec res( k_size );

  const double exp_scale = 1.0 / ( 2.0 * sigma * sigma );

  // Compute unnormalized kernel
  double sum = 0.0;
  for (size_t i = 0; i < k_size; ++i )
  {
    const double dx = ( static_cast<double>( i ) - static_cast<double>( half_k_size ) );
    res( i ) = exp( - dx * dx * exp_scale );
    sum += res( i );
  }

  // Normalize kernel
  const double inv = 1.0 / sum;
  for (size_t i = 0; i < k_size; ++i )
  {
    res( i ) *= inv;
  }

  return res;
}

/**
 ** Compute gaussian filtering of an image using user defined filter widths
 ** @param img Input image
 ** @param sigma standard deviation of kernel
 ** @param out Output image
 ** @param kernel_size_x Size of horizontal kernel (must be an odd number or 0 for automatic computation)
 ** @param kernel_size_y Size of vertical kernel (must be an add number or 0 for automatic computation)
 **/
template< typename Image >
void ImageGaussianFilter( const Image & img , const double sigma , Image & out ,
                          const size_t kernel_size_x , const size_t kernel_size_y )
{
  assert( kernel_size_x % 2 == 1 || kernel_size_x == 0 );
  assert( kernel_size_y % 2 == 1 || kernel_size_y == 0 );

  const Vec kernel_horiz = ComputeGaussianKernel( kernel_size_x , sigma );
  const Vec kernel_vert  = ComputeGaussianKernel( kernel_size_y , sigma );

  ImageSeparableConvolution( img , kernel_horiz , kernel_vert , out );
}

} // namespace image
} // namespace openMVG

#endif // OPENMVG_IMAGE_IMAGE_FILTERING_HPP
