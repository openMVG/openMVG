// Copyright (c) 2014 Romuald Perrot, Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_IMAGE_IMAGE_CONVOLUTION_CPP_OPENMP
#define OPENMVG_IMAGE_IMAGE_CONVOLUTION_CPP_OPENMP

#include <vector>
#include <Eigen/StdVector>

#include "openMVG/numeric/accumulator_trait.hpp"

namespace openMVG
{
  /**
   ** Horizontal convolution on a band of image
   ** - Assume kernel has odd size
   ** - Convolution on a band of an image
   ** - Assume out is already at correct size
   **
   ** @param img Input image
   ** @param kernel convolution kernel
   ** @param out Output image
   ** @param start_row Starting range of convolution range is [ start_row ; end_row [
   ** @param end_row End range of convolution [ start_row ; end_row [
   ** @param border_mgmt either BORDER_COPY or BORDER_CROP to tell what to do with borders
   **/
  template< typename Image, typename Kernel >
  void ImageHorizontalConvolution( const Image & img , const Kernel & kernel , Image & out ,
                                   const int start_row , const int end_row ,
                                   const EBorderManagement border_mgmt = BORDER_COPY )
  {
    typedef typename Image::Tpixel pix_t ;

    const int kernel_width = kernel.size() ;
    const int half_kernel_width = kernel_width / 2 ;

    const int cols = img.cols() ;

    if( border_mgmt == BORDER_COPY )
    {
      std::vector< pix_t, Eigen::aligned_allocator<pix_t> > line( cols + kernel_width ) ;

      for( int row = start_row ; row < end_row ; ++row )
      {
        // Copy line
        for( int k = 0 ; k < half_kernel_width ; ++k ) // pad before
        {
          line[ k ] = img.coeffRef( row , 0 ) ;
        }
        memcpy(&line[0] + half_kernel_width, img.data() + row * cols, sizeof(pix_t) * cols);
        for( int k = 0 ; k < half_kernel_width ; ++k ) // pad after
        {
          line[ k + half_kernel_width + cols ] = img.coeffRef( row , cols - 1 ) ;
        }

        // Apply convolution
        conv_buffer_( &line[0] , kernel.data(), cols , kernel_width );

        memcpy(out.data() + row * cols, &line[0], sizeof(pix_t) * cols);
      }
    }
    else
    {
      std::vector< pix_t, Eigen::aligned_allocator<pix_t> > line( cols ) ;

      for( int row = start_row ; row < end_row ; ++row )
      {
        // Copy line
        memcpy(&line[0], img.data() + row * cols, sizeof(pix_t) * cols);

        // Apply convolution
        for( int col = half_kernel_width ; col < cols - half_kernel_width ; ++col )
        {
          typename Accumulator<pix_t>::Type sum(0);
          for( int k = 0 ; k < kernel_width ; ++k )
          {
            sum += line[ k + col - half_kernel_width ] * kernel( k ) ;
          }
          out.coeffRef( row , col ) = sum ;
        }
      }
    }
  }


  /**
   ** Horizontal convolution on a band of image (in-place)
   ** assume kernel has odd size
   ** @param img Input image
   ** @param kernel convolution kernel
   ** @param out Output image
   ** @param start_row Starting range of convolution range is [ start_row ; end_row [
   ** @param end_row End range of convolution [ start_row ; end_row [
   ** @param border_mgmt either BORDER_COPY or BORDER_CROP to tell what to do with borders
   **/
  template< typename Image, typename Kernel >
  void ImageHorizontalConvolution( Image & img , const Kernel & kernel ,
                                   const int start_row , const int end_row ,
                                   const EBorderManagement border_mgmt = BORDER_COPY )
  {
    typedef typename Image::Tpixel pix_t ;

    const int kernel_width = kernel.size() ;
    const int half_kernel_width = kernel_width / 2 ;

    const int cols = img.cols() ;

    if( border_mgmt == BORDER_COPY )
    {
      std::vector< pix_t, Eigen::aligned_allocator<pix_t> > line( cols + kernel_width ) ;

      for( int row = start_row ; row < end_row ; ++row )
      {
        // Copy line
        for( int k = 0 ; k < half_kernel_width ; ++k ) // pad before
        {
          line[ k ] = img.coeffRef( row , 0 ) ;
        }
        memcpy(&line[0] + half_kernel_width, img.data() + row * cols, sizeof(pix_t) * cols);
        for( int k = 0 ; k < half_kernel_width ; ++k ) // pad after
        {
          line[ k + half_kernel_width + cols ] = img.coeffRef( row , cols - 1 ) ;
        }

        // Apply convolution
        conv_buffer_( &line[0], kernel.data(), cols , kernel_width );
        memcpy(img.data() + row * cols, &line[0], sizeof(pix_t) * cols);
      }
    }
    else
    {
      std::vector< pix_t, Eigen::aligned_allocator<pix_t> > line( cols ) ;

      for( int row = start_row ; row < end_row ; ++row )
      {
        // Copy line
        memcpy(&line[0], img.data() + row * cols, sizeof(pix_t) * cols);

        // Apply convolution
        for( int col = half_kernel_width ; col < cols - half_kernel_width ; ++col )
        {
          typename Accumulator<pix_t>::Type sum(0) ;
          for( int k = 0 ; k < kernel_width ; ++k )
          {
            sum += line[ k + col - half_kernel_width ] * kernel( k ) ;
          }
          img.coeffRef( row , col ) = sum ;
        }
      }
    }
  }


  /**
   ** Vertical (1d) convolution on a band (of column)
   ** assume out image is resized
   ** assume kernel has odd size
   ** @param img Input image
   ** @param kernel convolution kernel
   ** @param out Output image
   ** @param start_col Starting range of convolution range is [ start_col ; end_col [
   ** @param end_col End range of convolution [ start_col ; end_col [
   ** @param border_mgmt either BORDER_COPY or BORDER_CROP to tell what to do with borders
   **/
  template< typename Image, typename Kernel >
  void ImageVerticalConvolution( const Image & img , const Kernel & kernel , Image & out ,
                                 const int start_col , const int end_col ,
                                 const EBorderManagement border_mgmt = BORDER_COPY )
  {
    typedef typename Image::Tpixel pix_t ;

    const int kernel_width = kernel.size() ;
    const int half_kernel_width = kernel_width / 2 ;

    const int rows = img.rows() ;
    const int cols = img.cols() ;

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
        int col = start_col ;
        for( ; col < end_col - 4 ; col += 4 )
        {
          // Load Enough elements to compute convolution by kernel 
          // Compute first border

          sums[ 0 ] = sums[ 1 ] = sums[ 2 ] = sums[ 3 ] = typename Accumulator<pix_t>::Type(0) ; 

          for( int k = 0 ; k < kernel_width ; ++k )
          {
            const int id_row = input_index[ k ] ;
            const typename Accumulator<pix_t>::Type & kernel_value = kernel(k); 
            sums[ 0 ] += img.coeffRef( id_row , col ) * kernel_value ;
            sums[ 1 ] += img.coeffRef( id_row , col + 1 ) * kernel_value ;
            sums[ 2 ] += img.coeffRef( id_row , col + 2 ) * kernel_value ;
            sums[ 3 ] += img.coeffRef( id_row , col + 3 ) * kernel_value ;
          }

          out.coeffRef( row , col ) = sums[ 0 ] ;
          out.coeffRef( row , col + 1 ) = sums[ 1 ] ;
          out.coeffRef( row , col + 2 ) = sums[ 2 ] ; 
          out.coeffRef( row , col + 3 ) = sums[ 3 ] ;
        }

        // Compute last column elements    
        for( ; col < end_col ; ++col )
        {
          typename Accumulator<pix_t>::Type sum(0);

          for( int k = 0 ; k < kernel_width ; ++k )
          {
            const int id_row = input_index[ k ] ;
            const typename Accumulator<pix_t>::Type & kernel_value = kernel(k); 
            sum += img.coeffRef( id_row , col ) * kernel_value ;
          }

          out.coeffRef( row , col ) = sum ;
        }
      }
    }
    else
    {
      std::vector< pix_t, Eigen::aligned_allocator<pix_t>  > line( rows ) ;

      for( int col = start_col ; col < end_col ; ++col )
      {
        // Copy column
        for( int k = 0 ; k < rows ; ++k )
        {
          line[ k ] = img( k , col ) ;
        }

        // Apply convolution
        for( int row = half_kernel_width ; row < rows - half_kernel_width ; ++row )
        {
          typename Accumulator<pix_t>::Type sum(0);
          for( int k = 0 ; k < kernel_width ; ++k )
          {
            sum += line[ k + row - half_kernel_width ] * kernel( k ) ;
          }
          out.coeffRef( row , col ) = sum ;
        }
      }
    }
  }

  /**
   ** Vertical (1d) convolution of a band (in-place)
   ** assume kernel has odd size
   ** @param img Input/Output image
   ** @param kernel convolution kernel
   ** @param start_col Starting range of convolution range is [ start_col ; end_col [
   ** @param end_col End range of convolution [ start_col ; end_col [
   ** @param border_mgmt either BORDER_COPY or BORDER_CROP to tell what to do with borders
   **/
  template< typename Image, typename Kernel  >
  void ImageVerticalConvolution( Image & img , const Kernel & kernel ,
                                 const int start_col ,
                                 const int end_col ,
                                 const EBorderManagement border_mgmt = BORDER_COPY )
  {
    typedef typename Image::Tpixel pix_t ;

    const int kernel_width = kernel.size() ;
    const int half_kernel_width = kernel_width / 2 ;

    const int nb_row = img.rows() ;

    if( border_mgmt == BORDER_COPY )
    {
      std::vector< pix_t, Eigen::aligned_allocator<pix_t> > line( img.rows() + kernel_width ) ;

      for( int col = start_col ; col < end_col ; ++col )
      {
        // Copy column
        for( int k = 0 ; k < half_kernel_width ; ++k )
        {
          line[ k ] = img( 0 , col ) ;
        }
        for( int k = 0 ; k < nb_row ; ++k )
        {
          line[ k + half_kernel_width ] = img( k , col ) ;
        }
        for( int k = 0 ; k < half_kernel_width ; ++k )
        {
          line[ k + half_kernel_width + nb_row ] = img( nb_row - 1 , col ) ;
        }

        // Apply convolution
        conv_buffer_( &line[0], kernel.data() , nb_row , kernel_width );
        for( int row = 0 ; row < nb_row ; ++row )
        {
          img( row , col ) = line[row];
        }
      }
    }
    else
    {
      std::vector< pix_t, Eigen::aligned_allocator<pix_t> > line( img.rows() ) ;

      for( int col = start_col ; col < end_col ; ++col )
      {
        // Copy column
        for( int k = 0 ; k < img.rows() ; ++k )
        {
          line[ k ] = img( k , col ) ;
        }

        // Apply convolution
        for( int row = half_kernel_width ; row < nb_row - half_kernel_width ; ++row )
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
   ** Horizontal convolution using openMP threads
   ** assume kernel has odd size
   ** @param img Input image
   ** @param kernel convolution kernel
   ** @param out Output image
   ** @param border_mgmt either BORDER_COPY or BORDER_CROP to tell what to do with borders
   **/
  template< typename Image, typename Kernel >
  void ImageHorizontalConvolutionCPPThread( const Image & img , const Kernel & kernel , Image & out ,
      const EBorderManagement border_mgmt = BORDER_COPY )
  {
    // Resize element
    out.resize( img.Width() , img.Height() ) ;

    const int nb_thread = omp_get_max_threads();

    // Compute ranges
    std::vector< int > range;
    SplitRange( (int)0 , (int)img.rows() , nb_thread , range ) ;

#pragma omp parallel for
    for( int i = 1 ; i < static_cast<int>(range.size()) ; ++i )
      ImageHorizontalConvolution<Image, Kernel>(img, kernel, out, range[i-1], range[i], border_mgmt);
  }


  /**
   ** Horizontal convolution using openMP threads (in-place)
   ** assume kernel has odd size
   ** @param img Input/Output image
   ** @param kernel convolution kernel
   ** @param border_mgmt either BORDER_COPY or BORDER_CROP to tell what to do with borders
   **/
  template< typename Image, typename Kernel >
  void ImageHorizontalConvolutionCPPThread( Image & img , const Kernel & kernel ,
      const EBorderManagement border_mgmt = BORDER_COPY )
  {
    const int nb_thread = omp_get_max_threads();

    // Compute ranges
    std::vector< int > range;
    SplitRange( (int)0 , (int)img.rows() , nb_thread , range ) ;
    
#pragma omp parallel for
    for( int i = 1 ; i < static_cast<int>(range.size()) ; ++i )
      ImageHorizontalConvolution<Image, Kernel>(img, kernel, range[i-1], range[i], border_mgmt);
  }



  /**
   ** Vertical convolution using openMP threads
   ** assume kernel has odd size
   ** @param img Input image
   ** @param kernel convolution kernel
   ** @param out Output image
   ** @param border_mgmt either BORDER_COPY or BORDER_CROP to tell what to do with borders
   **/
  template< typename Image, typename Kernel >
  void ImageVerticalConvolutionCPPThread( const Image & img , const Kernel & kernel , Image & out ,
                                          const EBorderManagement border_mgmt = BORDER_COPY )
  {
    // Resize element
    out.resize( img.Width() , img.Height() ) ;

    const int nb_thread = omp_get_max_threads();

    // Compute ranges
    std::vector< int > range;
    SplitRange( (int)0 , (int)img.cols() , nb_thread , range ) ;

#pragma omp parallel for
    for( int i = 1 ; i < static_cast<int>(range.size()) ; ++i )
      ImageVerticalConvolution<Image, Kernel>(img, kernel, out, range[i-1], range[i], border_mgmt);
  }


  /**
   ** Vertical convolution using openMP threads (in-place)
   ** assume kernel has odd size
   ** @param img Input image
   ** @param kernel convolution kernel
   ** @param border_mgmt either BORDER_COPY or BORDER_CROP to tell what to do with borders
   **/
  template< typename Image, typename Kernel >
  void ImageVerticalConvolutionCPPThread( Image & img , const Kernel & kernel ,
                                          const EBorderManagement border_mgmt = BORDER_COPY )
  {
    const int nb_thread = omp_get_max_threads();

    // Compute ranges
    std::vector< int > range;
    SplitRange( (int)0 , (int)img.cols() , nb_thread , range ) ;

#pragma omp parallel for
    for( int i = 1 ; i < static_cast<int>(range.size()) ; ++i )
      ImageVerticalConvolution<Image, Kernel>(img, kernel, range[i-1], range[i], border_mgmt);
  }

} // namespace openMVG

#endif // OPENMVG_IMAGE_IMAGE_CONVOLUTION_CPP_OPENMP
