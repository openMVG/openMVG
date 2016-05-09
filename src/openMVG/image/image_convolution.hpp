// Copyright (c) 2014 Romuald Perrot, Pierre Moulon, Chris Sweeney.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_IMAGE_IMAGE_CONVOLUTION_HPP_
#define OPENMVG_IMAGE_IMAGE_CONVOLUTION_HPP_

#include "openMVG/numeric/numeric.h"
#include "openMVG/numeric/accumulator_trait.hpp"
#include "openMVG/image/image_container.hpp"

#include <cassert>
#include <vector>

/**
 ** @file Standard 2D image convolution functions :
 ** - vertical
 ** - horizontal
 ** - 2D (using standard 2d kernel of with separable kernels)
 **/

namespace openMVG
{
namespace image
{

/**
 ** General image convolution by a kernel
 ** assume kernel has odd size in both dimensions and (border pixel are copied)
 ** @param img source image
 ** @param kernel convolution kernel
 ** @param out resulting image
 **/
template< typename Image >
void ImageConvolution( const Image & img , const Mat & kernel , Image & out )
{
  const int kernel_width  = kernel.cols() ;
  const int kernel_height = kernel.rows() ;

  assert( kernel_width % 2 != 0 && kernel_height % 2 != 0 ) ;

  typedef typename Image::Tpixel pix_t ;
  typedef typename Accumulator< pix_t >::Type acc_pix_t ;

  out.resize( img.Width() , img.Height() ) ;

  for( int row = 0 ; row < img.rows() ; ++row )
  {
    for( int col = 0 ; col < img.cols() ; ++col )
    {
      acc_pix_t sum = acc_pix_t( ) ;

      for( int i = 0 ; i < kernel_height ; ++i )
      {
        for( int j = 0 ; j < kernel_width ; ++j )
        {
          int idy = row + i - kernel_height / 2 ;
          int idx = col + j - kernel_width / 2 ;

          // Get correct value
          idx = idx < 0 ? 0 : ( idx >= img.cols() ? img.cols() - 1 : idx ) ;
          idy = idy < 0 ? 0 : ( idy >= img.rows() ? img.rows() - 1 : idy ) ;

          sum += kernel( i , j ) * img( idy , idx ) ;
        }
      }
      out( row , col ) = sum ;
    }
  }
}

/**
 ** Horizontal (1d) convolution
 ** assume kernel has odd size
 ** @param img Input image
 ** @param kernel convolution kernel
 ** @param out Output image
 **/
template< typename ImageTypeIn , typename ImageTypeOut, typename Kernel >
void ImageHorizontalConvolution( const ImageTypeIn & img , const Kernel & kernel , ImageTypeOut & out )
{
  typedef typename ImageTypeIn::Tpixel pix_t ;

  const int rows ( img.rows() );
  const int cols ( img.cols() );

  out.resize( cols , rows ) ;

  const int kernel_width = kernel.size() ;
  const int half_kernel_width = kernel_width / 2 ;

  std::vector<pix_t, Eigen::aligned_allocator<pix_t> > line( cols + kernel_width );

  for( int row = 0 ; row < rows ; ++row )
  {
    // Copy line
    const pix_t start_pix = img.coeffRef( row , 0 ) ;
    for( int k = 0 ; k < half_kernel_width ; ++k ) // pad before
    {
      line[ k ] = start_pix ;
    }
    memcpy( &line[0] + half_kernel_width, img.data() + row * cols, sizeof( pix_t ) * cols );
    const pix_t end_pix = img.coeffRef( row , cols - 1 ) ;
    for( int k = 0 ; k < half_kernel_width ; ++k ) // pad after
    {
      line[ k + half_kernel_width + cols ] = end_pix ;
    }

    // Apply convolution
    conv_buffer_( &line[0] , kernel.data() , cols , kernel_width );

    memcpy( out.data() + row * cols, &line[0], sizeof( pix_t ) * cols );
  }
}

/**
 ** Vertical (1d) convolution
 ** assume kernel has odd size
 ** @param img Input image
 ** @param kernel convolution kernel
 ** @param out Output image
 **/
template< typename ImageTypeIn , typename ImageTypeOut, typename Kernel >
void ImageVerticalConvolution( const ImageTypeIn & img , const Kernel & kernel , ImageTypeOut & out )
{
  typedef typename ImageTypeIn::Tpixel pix_t ;

  const int kernel_width = kernel.size() ;
  const int half_kernel_width = kernel_width / 2 ;

  const int rows = img.rows() ;
  const int cols = img.cols() ;

  out.resize( cols , rows ) ;

  std::vector<pix_t, Eigen::aligned_allocator<pix_t> > line( rows + kernel_width );

  for( int col = 0 ; col < cols ; ++col )
  {
    // Copy column
    for( int k = 0 ; k < half_kernel_width ; ++k )
    {
      line[ k ] = img.coeffRef( 0 , col ) ;
    }
    for( int k = 0 ; k < rows ; ++k )
    {
      line[ k + half_kernel_width ] = img.coeffRef( k , col ) ;
    }
    for( int k = 0 ; k < half_kernel_width ; ++k )
    {
      line[ k + half_kernel_width + rows ] = img.coeffRef( rows - 1 , col ) ;
    }

    // Apply convolution
    conv_buffer_( &line[0] , kernel.data() , rows , kernel_width );

    for( int row = 0 ; row < rows ; ++row )
    {
      out.coeffRef( row , col ) = line[row];
    }
  }
}

/**
 ** Separable 2D convolution
 ** (nxm kernel is replaced by two 1D convolution of (size n then size m) )
 ** @param img source image
 ** @param horiz_k horizontal kernel
 ** @param vert_k vertical kernel
 ** @param out output image
 **/
template< typename ImageType, typename Kernel >
void ImageSeparableConvolution( const ImageType & img ,
                                const Kernel & horiz_k ,
                                const Kernel & vert_k ,
                                ImageType & out )
{
  // Cast the Kernel to the appropriate type
  typedef typename ImageType::Tpixel pix_t;
  typedef Eigen::Matrix<typename Accumulator<pix_t>::Type, Eigen::Dynamic, 1> VecKernel;
  const VecKernel horiz_k_cast = horiz_k.template cast< typename Accumulator<pix_t>::Type >();
  const VecKernel vert_k_cast = vert_k.template cast< typename Accumulator<pix_t>::Type >();

  ImageType tmp ;
  ImageHorizontalConvolution( img , horiz_k_cast , tmp ) ;
  ImageVerticalConvolution( tmp , vert_k_cast , out ) ;
}

typedef Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMatrixXf;

/// Specialization for Float based image (for arbitrary sized kernel)
inline void SeparableConvolution2d( const RowMatrixXf& image,
                                    const Eigen::Matrix<float, 1, Eigen::Dynamic>& kernel_x,
                                    const Eigen::Matrix<float, 1, Eigen::Dynamic>& kernel_y,
                                    RowMatrixXf* out )
{
  const int sigma_y = static_cast<int>( kernel_y.cols() );
  const int half_sigma_y = static_cast<int>( kernel_y.cols() ) / 2;

  // Convolving a vertical filter across rows is the same thing as transpose
  // multiply i.e. kernel_y^t * rows. This will give us the convoled value for
  // each row. However, care must be taken at the top and bottom borders.
  const Eigen::Matrix<float, 1, Eigen::Dynamic> reverse_kernel_y = kernel_y.reverse();
#if defined(OPENMVG_USE_OPENMP)
  #pragma omp parallel for schedule(dynamic)
#endif
  for ( int i = 0; i < half_sigma_y; i++ )
  {
    const int forward_size = i + half_sigma_y + 1;
    const int reverse_size = sigma_y - forward_size;
    out->row( i ) = kernel_y.tail( forward_size ) *
                    image.block( 0, 0, forward_size, image.cols() ) +
                    reverse_kernel_y.tail( reverse_size ) *
                    image.block( 1, 0, reverse_size, image.cols() );

    // Apply the same technique for the end rows.
    out->row( image.rows() - i - 1 ) =
      kernel_y.head( forward_size ) *
      image.block( image.rows() - forward_size, 0, forward_size, image.cols() )
      +
      reverse_kernel_y.head( reverse_size ) *
      image.block( image.rows() - reverse_size - 1, 0, reverse_size, image.cols() );
  }

  // Applying the rest of the y filter.
#if defined(OPENMVG_USE_OPENMP)
  #pragma omp parallel for schedule(dynamic)
#endif
  for ( int row = half_sigma_y; row < image.rows() - half_sigma_y; row++ )
  {
    out->row( row ) =  kernel_y * image.block( row - half_sigma_y, 0, sigma_y, out->cols() );
  }

  const int sigma_x = static_cast<int>( kernel_x.cols() );
  const int half_sigma_x = static_cast<int>( kernel_x.cols() / 2 );

  // Convolving with the horizontal filter is easy. Rather than using the kernel
  // as a sliding window, we use the row pixels as a sliding window around the
  // filter. We prepend and append the proper border values so that we are sure
  // to end up with the correct convolved values.
  Eigen::RowVectorXf temp_row( image.cols() + sigma_x - 1 );
#if defined(OPENMVG_USE_OPENMP)
  #pragma omp parallel for firstprivate(temp_row), schedule(dynamic)
#endif
  for ( int row = 0; row < out->rows(); row++ )
  {
    temp_row.head( half_sigma_x ) =
      out->row( row ).segment( 1, half_sigma_x ).reverse();
    temp_row.segment( half_sigma_x, image.cols() ) = out->row( row );
    temp_row.tail( half_sigma_x ) =
      out->row( row )
      .segment( image.cols() - 2 - half_sigma_x, half_sigma_x )
      .reverse();

    // Convolve the row. We perform the first step here explicitly so that we
    // avoid setting the row equal to zero.
    out->row( row ) = kernel_x( 0 ) * temp_row.head( image.cols() );
    for ( int i = 1; i < sigma_x; i++ )
    {
      out->row( row ) += kernel_x( i ) * temp_row.segment( i, image.cols() );
    }
  }
}


/**
* @brief Specialization for Image<float> in order to use SeparableConvolution2d
* @param img Input image
* @param horiz_k Kernel used for horizontal convolution
* @param vert_k Kernl used for vertical convolution
* @param[out] out Convolved image
*/
template<typename Kernel>
void ImageSeparableConvolution( const Image<float> & img ,
                                const Kernel & horiz_k ,
                                const Kernel & vert_k ,
                                Image<float> & out )
{
  // Cast the Kernel to the appropriate type
  typedef Image<float>::Tpixel pix_t;
  typedef Eigen::Matrix<typename openMVG::Accumulator<pix_t>::Type, Eigen::Dynamic, 1> VecKernel;
  const VecKernel horiz_k_cast = horiz_k.template cast< typename openMVG::Accumulator<pix_t>::Type >();
  const VecKernel vert_k_cast = vert_k.template cast< typename openMVG::Accumulator<pix_t>::Type >();

  out.resize( img.Width(), img.Height() );
  SeparableConvolution2d( img.GetMat(), horiz_k_cast, vert_k_cast, &( ( Image<float>::Base& )out ) );
}

} // namespace image
} // namespace openMVG

#endif // OPENMVG_IMAGE_IMAGE_CONVOLUTION_HPP_
