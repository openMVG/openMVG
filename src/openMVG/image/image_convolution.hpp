// Copyright (c) 2014 Romuald Perrot, Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_IMAGE_IMAGE_CONVOLUTION_HPP_
#define OPENMVG_IMAGE_IMAGE_CONVOLUTION_HPP_

#include "openMVG/numeric/accumulator_trait.hpp"
/**
 ** @file Standard 2D image convolution functions :
 ** - vertical (with and without in place modifications)
 ** - horizontal (with and without in place modifications)
 ** - 2D (using standard 2d kernel of with separable kernels -- both with and without in place modifications)
 **/

namespace openMVG
{
  /**
   ** General image convolution by a kernel
   ** assume kernel has odd size in both dimensions
   ** @param img source image
   ** @param kernel convolution kernel
   ** @param out resulting image
   ** @param border_mgmt either BORDER_COPY or BORDER_CROP to tell what to do with borders
   **/
  template< typename Image >
  void ImageConvolution( const Image & img , const Mat & kernel , Image & out , const EBorderManagement border_mgmt = BORDER_COPY )
  {
    const int kernel_width  = kernel.cols() ;
    const int kernel_height = kernel.rows() ;

    assert( kernel_width % 2 != 0 && kernel_height % 2 != 0 ) ;

    typedef typename Image::Tpixel pix_t ;
    typedef typename Accumulator< pix_t >::Type acc_pix_t ;

    out.resize( img.Width() , img.Height() ) ;
    if( border_mgmt == BORDER_COPY )
    {
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
    else // border CROP
    {
      for( int row = kernel_height / 2 ; row < img.rows() - kernel_height / 2 ; ++row )
      {
        for( int col = kernel_width / 2 ; col < img.cols() - kernel_width / 2 ; ++col )
        {
          acc_pix_t sum(0) ;
          // Compute value for
          for( int i = 0 ; i < kernel_height ; ++i )
          {
            for( int j = 0 ; j < kernel_width ; ++j )
            {
              sum += kernel( i , j ) * img( row + i - kernel_height / 2 , col + j - kernel_width / 2 ) ;
            }
          }
          out( row , col ) = sum ;
        }
      }
    }
  }


  /**
   ** General image convolution by a kernel (in place)
   ** assume kernel has odd size in both dimensions
   ** @param img image source and destination image
   ** @param kernel convolution kernel
   ** @param border_mgmt either BORDER_COPY or BORDER_CROP to tell what to do with borders
   **/
  template< typename Image >
  void ImageConvolution( Image & img , const Mat & kernel , const EBorderManagement border_mgmt = BORDER_COPY )
  {
    const int kernel_width  = kernel.cols() ;
    const int kernel_height = kernel.rows() ;

    assert( kernel_width % 2 != 0 && kernel_height % 2 != 0 ) ;

    typedef typename Image::Tpixel pix_t ;
    typedef typename Accumulator< pix_t >::Type acc_pix_t ;

    Image out ;
    out.resize( img.Width() , img.Height() ) ;
    if( border_mgmt == BORDER_COPY )
    {
      for( int row = 0 ; row < img.rows() ; ++row )
      {
        for( int col = 0 ; col < img.cols() ; ++col )
        {
          pix_t sum = pix_t( ) ;

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
    else // border CROP
    {
      for( int row = kernel_height / 2 ; row < img.rows() - kernel_height / 2 ; ++row )
      {
        for( int col = kernel_width / 2 ; col < img.cols() - kernel_width / 2 ; ++col )
        {
          acc_pix_t sum(0);
          // Compute value for
          for( int i = 0 ; i < kernel_height ; ++i )
          {
            for( int j = 0 ; j < kernel_width ; ++j )
            {
              sum += kernel( i , j ) * img( row + i - kernel_height / 2 , col + j - kernel_width / 2 ) ;
            }
          }
          out( row , col ) = sum ;
        }
      }
    }
    img = out ;
  }


  /**
   ** Horizontal (1d) convolution
   ** assume kernel has odd size
   ** @param img Input image
   ** @param kernel convolution kernel
   ** @param out Output image
   ** @param border_mgmt either BORDER_COPY or BORDER_CROP to tell what to do with borders
   **/
  template< typename ImageTypeIn , typename ImageTypeOut, typename Kernel >
  static void ImageHorizontalConvolution( const ImageTypeIn & img , const Kernel & kernel , ImageTypeOut & out , const EBorderManagement border_mgmt = BORDER_COPY)
  {
    typedef typename ImageTypeIn::Tpixel pix_t ;

    const int rows ( img.rows() );
    const int cols ( img.cols() );

    out.resize( cols , rows ) ;

    const int kernel_width = kernel.size() ;
    const int half_kernel_width = kernel_width / 2 ;

    if( border_mgmt == BORDER_COPY )
    {
      std::vector<pix_t, Eigen::aligned_allocator<pix_t> > line( cols + kernel_width );

      for( int row = 0 ; row < rows ; ++row )
      {
        // Copy line
        const pix_t start_pix = img.coeffRef( row , 0 ) ;
        for( int k = 0 ; k < half_kernel_width ; ++k ) // pad before
        {
          line[ k ] = start_pix ;
        }
        memcpy(&line[0] + half_kernel_width, img.data() + row * cols, sizeof(pix_t) * cols);
        const pix_t end_pix = img.coeffRef( row , cols - 1 ) ;
        for( int k = 0 ; k < half_kernel_width ; ++k ) // pad after
        {
          line[ k + half_kernel_width + cols ] = end_pix ;
        }

        // Apply convolution
        conv_buffer_( &line[0] , kernel.data() , cols , kernel_width );

        memcpy(out.data() + row * cols, &line[0], sizeof(pix_t) * cols);
      }
    }
    else
    {
      std::vector< pix_t, Eigen::aligned_allocator<pix_t>  > line( cols ) ;
      typedef typename Accumulator< pix_t >::Type acc_pix_t ;
      
      for( int row = 0 ; row < rows ; ++row )
      {
        // Copy line
        memcpy(&line[0], img.data() + row * cols, sizeof(pix_t) * cols);

        // Apply convolution
        for( int col = half_kernel_width ; col < cols - half_kernel_width ; ++col )
        {
          acc_pix_t sum(0);
          for( int k = 0 ; k < kernel_width ; ++k )
          {
            sum += line[ k + col - half_kernel_width ] * static_cast<pix_t>( kernel( k ) );
          }
          out( row , col ) = sum ;
        }
      }
    }
  }

  /**
   ** Horizontal (1d) convolution (in-place)
   ** assume kernel has odd size
   ** @param img Input image
   ** @param kernel convolution kernel
   ** @param out Output image
   ** @param border_mgmt either BORDER_COPY or BORDER_CROP to tell what to do with borders
   **/
  template< typename Image, typename Kernel >
  static void ImageHorizontalConvolution( Image & img , const Kernel & kernel , const EBorderManagement border_mgmt = BORDER_COPY )
  {
    typedef typename Image::Tpixel pix_t ;

    const int rows = img.rows() ;
    const int cols = img.cols() ;

    const int kernel_width = kernel.size() ;
    const int half_kernel_width = kernel_width / 2 ;

    if( border_mgmt == BORDER_COPY )
    {
      std::vector< pix_t, Eigen::aligned_allocator<pix_t>  > line( cols + kernel_width ) ;

      for( int row = 0 ; row < rows ; ++row )
      {
        // Copy line
        for( int k = 0 ; k < half_kernel_width ; ++k ) // pad before
        {
          line[ k ] = img.coeffRef( row , 0 ) ;
        }
        for( int k = 0 ; k < cols ; ++k )
        {
          line[ k + half_kernel_width ] = img.coeffRef( row , k ) ;
        }
        for( int k = 0 ; k < half_kernel_width ; ++k ) // pad after
        {
          line[ k + half_kernel_width + cols ] = img.coeffRef( row , cols - 1 ) ;
        }

        // Apply convolution
        conv_buffer_( &line[0], kernel.data(), cols , kernel_width );

        for( int col = 0 ; col < cols ; ++col )
        {
          img( row , col ) = line[col];
        }
      }
    }
    else
    {
      std::vector< pix_t, Eigen::aligned_allocator<pix_t>  > line( img.cols() ) ;

      for( int row = half_kernel_width ; row < img.rows() - half_kernel_width ; ++row )
      {
        // Copy line
        for( int k = 0 ; k < img.cols() ; ++k )
        {
          line[ k ] = img( row , k ) ;
        }

        // Apply convolution
        for( int col = half_kernel_width ; col < img.cols() - half_kernel_width ; ++col )
        {
          typename Accumulator<pix_t>::Type sum(0) ;
          for( int k = 0 ; k < kernel_width ; ++k )
          {
            sum += line[ k + col - half_kernel_width ] * kernel( k ) ;
          }
          img( row , col ) = sum ;
        }
      }
    }
  }


  /**
   ** Vertical (1d) convolution
   ** assume kernel has odd size
   ** @param img Input image
   ** @param kernel convolution kernel
   ** @param out Output image
   ** @param border_mgmt either BORDER_COPY or BORDER_CROP to tell what to do with borders
   **/
  template< typename ImageTypeIn , typename ImageTypeOut, typename Kernel >
  void ImageVerticalConvolution( const ImageTypeIn & img , const Kernel & kernel , ImageTypeOut & out , const EBorderManagement border_mgmt = BORDER_COPY)
  {
    typedef typename ImageTypeIn::Tpixel pix_t ;

    const int kernel_width = kernel.size() ;
    const int half_kernel_width = kernel_width / 2 ;

    const int rows = img.rows() ;
    const int cols = img.cols() ;

    out.resize( cols , rows ) ;

    if( border_mgmt == BORDER_COPY )
    {
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
    else
    {
      std::vector< pix_t, Eigen::aligned_allocator<pix_t>  > line( img.rows() ) ;

      for( int col = 0 ; col < img.cols() ; ++col )
      {
        // Copy column
        for( int k = 0 ; k < img.rows() ; ++k )
        {
          line[ k ] = img( k , col ) ;
        }

        // Apply convolution
        for( int row = half_kernel_width ; row < img.rows() - half_kernel_width ; ++row )
        {
          typename Accumulator<pix_t>::Type sum(0) ;
          for( int k = 0 ; k < kernel_width ; ++k )
          {
            sum += line[ k + row - half_kernel_width ] * kernel( k ) ;
          }
          out( row , col ) = sum ;
        }
      }
    }
  }

  /**
   ** Vertical (1d) convolution (in-place)
   ** assume kernel has odd size
   ** @param img Input/Output image
   ** @param kernel convolution kernel
   ** @param border_mgmt either BORDER_COPY or BORDER_CROP to tell what to do with borders
   **/
  template< typename Image, typename Kernel >
  void ImageVerticalConvolution( Image & img , const Kernel & kernel , const EBorderManagement border_mgmt = BORDER_COPY )
  {
    typedef typename Image::Tpixel pix_t ;

    const int kernel_width = kernel.size() ;
    const int half_kernel_width = kernel_width / 2 ;

    if( border_mgmt == BORDER_COPY )
    {
      std::vector< pix_t, Eigen::aligned_allocator<pix_t> > line ;
      line.resize( img.rows() + kernel_width ) ;

      for( int col = 0 ; col < img.cols() ; ++col )
      {
        // Copy column
        for( int k = 0 ; k < half_kernel_width ; ++k )
        {
          line[ k ] = img( 0 , col ) ;
        }
        for( int k = 0 ; k < img.rows() ; ++k )
        {
          line[ k + half_kernel_width ] = img( k , col ) ;
        }
        for( int k = 0 ; k < half_kernel_width ; ++k )
        {
          line[ k + half_kernel_width + img.rows() ] = img( img.rows() - 1 , col ) ;
        }

        // Apply convolution
        conv_buffer_( &line[0], kernel.data() , img.rows(), kernel_width );
        for( int row = 0 ; row < img.rows() ; ++row )
        {
          img( row , col ) = line[row];
        }
      }
    }
    else
    {
      std::vector< pix_t, Eigen::aligned_allocator<pix_t> > line( img.rows() ) ;

      for( int col = 0 ; col < img.cols() ; ++col )
      {
        // Copy column
        for( int k = 0 ; k < img.rows() ; ++k )
        {
          line[ k ] = img( k , col ) ;
        }

        // Apply convolution
        for( int row = half_kernel_width ; row < img.rows() - half_kernel_width ; ++row )
        {
          typename Accumulator<pix_t>::Type sum(0) ;
          for( int k = 0 ; k < kernel_width ; ++k )
          {
            sum += line[ k + row - half_kernel_width ] * kernel( k ) ;
          }
          img( row , col ) = sum ;
        }
      }
    }
  }

  /**
   ** Vertical (1d) convolution (in-place)
   ** assume kernel has odd size
   ** @param img Input/Output image
   ** @param kernel convolution kernel
   ** @param border_mgmt either BORDER_COPY or BORDER_CROP to tell what to do with borders
   **/
  template< typename Image, typename Kernel >
  void ImageVerticalConvolution( const Image & img , const Kernel & kernel , Image & out , const EBorderManagement border_mgmt = BORDER_COPY )
  {
    typedef typename Image::Tpixel pix_t ;

    const int kernel_width = kernel.size() ;
    const int half_kernel_width = kernel_width / 2 ;

    const int rows ( img.rows() );
    const int cols ( img.cols() );

    out.resize( cols , rows ) ;

    if( border_mgmt == BORDER_COPY )
    {
      typename Accumulator<pix_t>::Type sums[ 4 ] ;

      std::vector<int> input_index( kernel_width );

      for( int row = 0 ; row < rows ; ++row )
      {
        // Compute prev and next index (in row) around pixel
        for( int i = 0 ; i < kernel_width ; ++i )
        {
          input_index[ i ] = ( row - half_kernel_width + i ) ;
          input_index[ i ] = std::max( 0 , input_index[ i ] ) ;
          input_index[ i ] = std::min( rows - 1 , input_index[ i ] ) ;
        }

        // Compute using 4 contiguous values
        int col = 0 ;
        for( ; col < cols - 4 ; col += 4 )
        {
          // Load Enough elements to compute convolution by kernel
          // Compute first border

          sums[ 0 ] = sums[ 1 ] = sums[ 2 ] = sums[ 3 ] = typename Accumulator<pix_t>::Type(0) ;

          for( int k = 0 ; k < kernel_width ; ++k )
          {
            const int id_row = input_index[ k ] ;
            sums[ 0 ] += img.coeffRef( id_row , col ) * kernel(k) ;
            sums[ 1 ] += img.coeffRef( id_row , col + 1 ) * kernel(k) ;
            sums[ 2 ] += img.coeffRef( id_row , col + 2 ) * kernel(k) ;
            sums[ 3 ] += img.coeffRef( id_row , col + 3 ) * kernel(k) ;
          }

          out.coeffRef( row , col ) = sums[ 0 ] ;
          out.coeffRef( row , col + 1 ) = sums[ 1 ] ;
          out.coeffRef( row , col + 2 ) = sums[ 2 ] ;
          out.coeffRef( row , col + 3 ) = sums[ 3 ] ;
        }

        // Compute last column elements
        for( ; col < cols ; ++col )
        {
          typename Accumulator<pix_t>::Type sum(0) ;

          for( int k = 0 ; k < kernel_width ; ++k )
          {
            const int id_row = input_index[ k ] ;
            const typename Accumulator<pix_t>::Type & kernel_value = kernel(k) ;
            sum += img.coeffRef( id_row , col ) * kernel_value ;
          }
          out.coeffRef( row , col ) = sum ;
        }
      }
    }
    else
    {
      std::vector< pix_t, Eigen::aligned_allocator<pix_t>  > line( rows ) ;

      for( int col = 0 ; col < cols ; ++col )
      {
        // Copy column
        for( int k = 0 ; k < rows ; ++k )
        {
          line[ k ] = img( k , col ) ;
        }

        // Apply convolution
        for( int row = half_kernel_width ; row < rows - half_kernel_width ; ++row )
        {
          typename Accumulator<pix_t>::Type sum(0) ;
          for( int k = 0 ; k < kernel_width ; ++k )
          {
            sum += line[ k + row - half_kernel_width ] * static_cast<pix_t>( kernel( k ) );
          }
          out( row , col ) = sum ;
        }
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
   ** @param border_mgmt either BORDER_COPY or BORDER_CROP to tell what to do with borders
   **/
  template< typename ImageType, typename Kernel >
  void ImageSeparableConvolution( const ImageType & img ,
                                  const Kernel & horiz_k ,
                                  const Kernel & vert_k ,
                                  ImageType & out ,
                                  const EBorderManagement border_mgmt = BORDER_COPY )
  {
    // Cast the Kernel to the appropriate type
    typedef typename ImageType::Tpixel pix_t;
    typedef Eigen::Matrix<typename Accumulator<pix_t>::Type, Eigen::Dynamic, 1> VecKernel;
    const VecKernel horiz_k_cast = horiz_k.template cast< typename Accumulator<pix_t>::Type >();
    const VecKernel vert_k_cast = vert_k.template cast< typename Accumulator<pix_t>::Type >();
    
    ImageType tmp ;
#if defined(USE_CPP11) && defined(HAVE_CXX11_THREAD)
    ImageHorizontalConvolutionCPPThread( img , horiz_k_cast , tmp , border_mgmt ) ;
    ImageVerticalConvolutionCPPThread( tmp , vert_k_cast , out , border_mgmt ) ;
#else
#if defined(USE_OPENMP)
    ImageHorizontalConvolutionCPPThread( img , horiz_k_cast , tmp , border_mgmt ) ;
    ImageVerticalConvolutionCPPThread( tmp , vert_k_cast , out , border_mgmt ) ;
#else
    ImageHorizontalConvolution( img , horiz_k_cast , tmp , border_mgmt ) ;
    ImageVerticalConvolution( tmp , vert_k_cast , out , border_mgmt ) ;
#endif
#endif
  }


  /**
   ** Separable 2D convolution (in place)
   ** (nxm kernel is replaced by two 1D convolution of (size n then size m) )
   ** @param img image source and destination image
   ** @param horiz_k horizontal kernel
   ** @param vert_k vertical kernel
   ** @param border_mgmt either BORDER_COPY or BORDER_CROP to tell what to do with borders
   **/
  template< typename Image, typename Kernel >
  void ImageSeparableConvolution( Image & img ,
                                  const Kernel & horiz_k ,
                                  const Kernel & vert_k ,
                                  const EBorderManagement border_mgmt = BORDER_COPY )
  {
    // Cast the Kernel to the appropriate type
    typedef typename Image::Tpixel pix_t;
    typedef Eigen::Matrix<typename Accumulator<pix_t>::Type, Eigen::Dynamic, 1> VecKernel;
    const VecKernel horiz_k_cast = horiz_k.template cast< typename Accumulator<pix_t>::Type >();
    const VecKernel vert_k_cast = vert_k.template cast< typename Accumulator<pix_t>::Type >();

#if defined(USE_CPP11) && defined(HAVE_CXX11_THREAD)
    ImageHorizontalConvolutionCPPThread( img , horiz_k_cast , border_mgmt ) ;
    ImageVerticalConvolutionCPPThread( img , vert_k_cast , border_mgmt ) ;
#else
#if defined(USE_OPENMP)
    ImageHorizontalConvolutionCPPThread( img , horiz_k_cast , border_mgmt ) ;
    ImageVerticalConvolutionCPPThread( img , vert_k_cast , border_mgmt ) ;
#else
    ImageHorizontalConvolution( img , horiz_k_cast , border_mgmt ) ;
    ImageVerticalConvolution( img , vert_k_cast , border_mgmt ) ;
#endif
#endif
  }


} // namespace openMVG

#endif // OPENMVG_IMAGE_IMAGE_CONVOLUTION_HPP_

