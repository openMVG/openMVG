// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_IMAGE_IMAGE_IMAGE_IO_HPP
#define OPENMVG_IMAGE_IMAGE_IMAGE_IO_HPP

#include "openMVG/image/image_container.hpp"
#include "openMVG/image/image_converter.hpp"

#include <cstdio>
#include <vector>

namespace openMVG
{
namespace image
{

extern template class Image<unsigned char>;
extern template class Image<float>;
extern template class Image<double>;
extern template class Image<RGBColor>;
extern template class Image<RGBAColor>;

/**
* @enum Format List of input format handled by IO operations
* @var Pnm
*   Generic PBM, PGM or PPM format
* @var Png
*   Portable Network Graphics
* @var Jpg
*   Joint Photographic Experts Group
* @var Tiff
*   Tagged Image File Format
*/
enum Format
{
  Pnm, Png, Jpg, Tiff, Unknown
};


/**
* @brief Get format of an image
* @param c Input path of the image
* @return Image format of the specified image
*/
Format GetFormat( const char *c );



/**
* @brief Load an image<T> from the provided input filename
* @param path Input path of the image to load
* @param[out] Output image
* @retval 1 If loading is correct
* @retval 0 If there was an error during load operation
*/
template<typename T>
int ReadImage( const char * path , Image<T> * image );

/**
* @brief Save an image<T> from the provided input filename
* @param path Output path of the image to save
* @param image Image to save
* @retval 1 If saving is correct
* @retval 0 If there was an error during save operation
*/
template<typename T>
int WriteImage( const char * path , const Image<T>& image );

/**
* @brief Unsigned char specialization (The memory pointer must be null as input)
* @param path Input path of the image to load
* @param[out] image Output image
* @param[out] w Width of the loaded image
* @param[out] h Height of the loaded image
* @param[out] depth Depth of the image
* @retval 1 If loading is correct
* @retval 0 If there was an error during load operation
*/
int ReadImage( const char * path, std::vector<unsigned char> * image , int * w, int * h, int * depth );

/**
* @brief Unsigned char specialization
* @param path Output path of the image to save
* @param array Image to save
* @param w Width of the image to save
* @param h Height of the image to save
* @param depth Depth of the image to save
* @retval 1 If saving is correct
* @retval 0 If there was an error during save operation
*/
int WriteImage( const char * path , const std::vector<unsigned char>& array, const int w, const int h, const int depth );

//--
// PNG I/O
//--
/**
* @brief Read PNG file from file
* @param[in] path Input file path
* @param[out] array Output image data array
* @param[out] w Image width
* @param[out] h Image height
* @param[out] depth Depth of image
* @retval 0 if there was an error during read operation
* @return non nul value if read operation is valid
*/
int ReadPng( const char * path , std::vector<unsigned char> * array , int * w, int * h, int * depth );

/**
* @brief Read PNG file from a stream
* @param[in] stream Input data stream
* @param[out] array Output image data array
* @param[out] w Image width
* @param[out] h Image height
* @param[out] depth Depth of image
* @retval 0 if there was an error during read operation
* @return non nul value if read operation is valid
*/
int ReadPngStream( FILE * stream , std::vector<unsigned char> * array , int * w, int * h, int * depth );


/**
* @brief Write PNG file to a file
* @param path Output file path
* @param array input image data array
* @param w Image width
* @param h Image height
* @param depth Depth of image
* @retval 0 if there was an error during write operation
* @return non nul value if write operation is valid
*/
int WritePng( const char * path , const std::vector<unsigned char>& array, int w, int h, int depth );

/**
* @brief Write PNG file to a stream
* @param[out] stream Output data stream
* @param array input image data array
* @param w Image width
* @param h Image height
* @param depth Depth of image
* @retval 0 if there was an error during write operation
* @return non nul value if write operation is valid
*/
int WritePngStream( FILE * stream ,  const std::vector<unsigned char>& array, int w, int h, int depth );

//--
// JPG I/O
//--
/**
* @brief Read JPEG image from file
* @param[in] path Input filepath
* @param[out] array Output image data
* @param[out] w Image width
* @param[out] h Image height
* @param[out] depth Depth of image
* @retval 0 if there is an error during read operation
* @return non nul value if read operation is valid
*/
int ReadJpg( const char * path , std::vector<unsigned char> * array, int * w, int * h, int * depth );

/**
* @brief Read JPEG image from stream
* @param[in] stream Input data stream
* @param[out] array Output image data
* @param[out] w Image width
* @param[out] h Image height
* @param[out] depth Depth of image
* @retval 0 if there is an error during read operation
* @return non nul value if read operation is valid
*/
int ReadJpgStream( FILE * stream , std::vector<unsigned char> * array, int * w, int * h, int * depth );

/**
* @brief Write JPEG file
* @param path Output image path
* @param img Input image
* @param quality JPG quality
* @retval 0 if there is an error during write operation
* @return non nul value if write operation is valid
* @note Quality is between 0 to 100
*/
template<typename T>
int WriteJpg( const char * path , const Image<T>& img , int quality = 90 );

/**
* @brief Write JPEG file
* @param path Output image path
* @param array Input image data
* @param w Width of the image
* @param h Height of the image
* @param depth Depth of the image
* @param quality JPG quality
* @retval 0 if there is an error during write operation
* @return non nul value if write operation is valid
* @note Quality is between 0 to 100
*/
int WriteJpg( const char * path , const std::vector<unsigned char>& array, int w, int h, int depth, int quality = 90 );

/**
* @brief Write JPEG file to a stream
* @param[out] stream Output data stream
* @param array Input image data
* @param w Width of the image
* @param h Height of the image
* @param depth Depth of the image
* @param quality JPG quality
* @retval 0 if there is an error during write operation
* @return non nul value if write operation is valid
* @note Quality is between 0 to 100
*/
int WriteJpgStream( FILE * stream , const std::vector<unsigned char>& array, int w, int h, int depth, int quality = 90 );


//--
// PNM/PGM I/O
//--
/**
* @brief Read PNM/PGM file
* @param[in] path Input image path
* @param[out] array Output image data
* @param[out] w Width of the image
* @param[out] h Height of the image
* @param[out] depth Depth of the image
* @retval 0 if there was an error during read operation
* @return non nul value if there was an error during read operation
*/
int ReadPnm( const char * path , std::vector<unsigned char> * array, int * w, int * h, int * depth );

/**
* @brief Read PNM/PGM from a stream
* @param[in] stream Input image stream data
* @param[out] array Output image data
* @param[out] w Width of the image
* @param[out] h Height of the image
* @param[out] depth Depth of the image
* @retval 0 if there was an error during read operation
* @return non nul value if there was an error during read operation
*/
int ReadPnmStream( FILE * stream , std::vector<unsigned char> * array, int * w, int * h, int * depth );

/**
* @brief Write PNM/PGM from to a file
* @param[in] path Output image path
* @param array input image data
* @param w Width of the image
* @param h Height of the image
* @param depth Depth of the image
* @retval 0 if there was an error during write operation
* @return non nul value if there was an error during write operation
*/
int WritePnm( const char * path , const std::vector<unsigned char>& array, int w, int h, int depth );

/**
* @brief Write PNM/PGM to a stream
* @param[in] stream Output image stream data
* @param array input image data
* @param w Width of the image
* @param h Height of the image
* @param depth Depth of the image
* @retval 0 if there was an error during write operation
* @return non nul value if there was an error during write operation
*/
int WritePnmStream( FILE *,  const std::vector<unsigned char>& array, int w, int h, int depth );

//--
// TIFF I/O
//--
/**
* @brief Read TIFF image from a file
* @param path Input file path
* @param[out] array Output image data
* @param[out] w Width of the image
* @param[out] h Height of the image
* @param[out] depth Depth of the image
* @retval 0 if there was an error during read operation
* @return non nul value if there was an error during read operation
*/
int ReadTiff( const char * path , std::vector<unsigned char> * array, int * w, int * h, int * depth );

/**
* @brief write TIFF image to a file
* @param path Output file path
* @param array Output image data
* @param w Width of the image
* @param h Height of the image
* @param depth Depth of the image
* @retval 0 if there was an error during write operation
* @return non nul value if there was an error during write operation
*/
int WriteTiff( const char * path , const std::vector<unsigned char>& array, int w, int h, int depth );

/**
* @brief  structure used to know the size of an image
* @note: can be extented later to support the pixel data type and the number of channel:
* i.e: - unsigned char, 1 => gray image
*      - unsigned char, 3 => rgb image
*      - unsigned char, 4 => rgba image
*/
struct ImageHeader
{
  /// Width of the image
  int width;

  /// Height of the image
  int height;
};

/**
* @brief Read image header from a file
* @param path Input image file path
* @param[out] hdr Output header
* @retval true If header is correctly read
* @retval false If header could not be read
*/
bool ReadImageHeader( const char * path , ImageHeader * hdr );

/**
* @brief Read PNG image header from a file
* @param path Input image file path
* @param[out] hdr Output header
* @retval true If header is correctly read
* @retval false If header could not be read
*/
bool Read_PNG_ImageHeader( const char * path , ImageHeader * hdr );

/**
* @brief Read JPEG image header from a file
* @param path Input image file path
* @param[out] hdr Output header
* @retval true If header is correctly read
* @retval false If header could not be read
*/
bool Read_JPG_ImageHeader( const char * path , ImageHeader * hdr );

/**
* @brief Read PNM/PGM image header from a file
* @param path Input image file path
* @param[out] hdr Output header
* @retval true If header is correctly read
* @retval false If header could not be read
*/
bool Read_PNM_ImageHeader( const char * path , ImageHeader * hdr );

/**
* @brief Read TIFF image header from a file
* @param path Input image file path
* @param[out] hdr Output header
* @retval true If header is correctly read
* @retval false If header could not be read
*/
bool Read_TIFF_ImageHeader( const char * path , ImageHeader * hdr );


/**
* @brief Generic Image read from file
* @param[in] path Input image path
* @param[out] im Ouput image
* @retval 0 if there was an errir during read operation
* @retval 1 if read is correct
*/
template<>
inline int ReadImage( const char * path, Image<unsigned char> * im )
{
  std::vector<unsigned char> ptr;
  int w, h, depth;
  const int res = ReadImage( path, &ptr, &w, &h, &depth );
  if ( res == 1 && depth == 1 )
  {
    //convert raw array to Image
    ( *im ) = Eigen::Map<Image<unsigned char>::Base>( &ptr[0], h, w );
  }
  else if ( res == 1 && depth == 3 )
  {
    //-- Must convert RGB to gray
    RGBColor * ptrCol = reinterpret_cast<RGBColor*>( &ptr[0] );
    Image<RGBColor> rgbColIm;
    rgbColIm = Eigen::Map<Image<RGBColor>::Base>( ptrCol, h, w );
    //convert RGB to gray
    ConvertPixelType( rgbColIm, im );
  }
  else if ( res == 1 && depth == 4 )
  {
    //-- Must convert RGBA to gray
    RGBAColor * ptrCol = reinterpret_cast<RGBAColor*>( &ptr[0] );
    Image<RGBAColor> rgbaColIm;
    rgbaColIm = Eigen::Map<Image<RGBAColor>::Base>( ptrCol, h, w );
    //convert RGBA to gray
    ConvertPixelType( rgbaColIm, im );
  }
  else if ( depth != 1 )
  {
    return 0;
  }
  return res;
}


/**
* @brief Generic Image read from file (overload for RGBColor)
* @param[in] path Input image path
* @param[out] im Ouput image
* @retval 0 if there was an errir during read operation
* @retval 1 if read is correct
*/
template<>
inline int ReadImage( const char * path, Image<RGBColor> * im )
{
  std::vector<unsigned char> ptr;
  int w, h, depth;
  const int res = ReadImage( path, &ptr, &w, &h, &depth );
  if ( res == 1 && depth == 3 )
  {
    RGBColor * ptrCol = reinterpret_cast<RGBColor*>( &ptr[0] );
    //convert raw array to Image
    ( *im ) = Eigen::Map<Image<RGBColor>::Base>( ptrCol, h, w );
  }
  else if ( res == 1 && depth == 4 )
  {
    //-- Must convert RGBA to RGB
    RGBAColor * ptrCol = reinterpret_cast<RGBAColor*>( &ptr[0] );
    Image<RGBAColor> rgbaColIm;
    rgbaColIm = Eigen::Map<Image<RGBAColor>::Base>( ptrCol, h, w );
    //convert RGBA to RGB
    ConvertPixelType( rgbaColIm, im );
  }
  else
  {
    return 0;
  }
  return res;
}

/**
* @brief Generic Image read from file (overload for RGBAColor)
* @param[in] path Input image path
* @param[out] im Ouput image
* @retval 0 if there was an errir during read operation
* @retval 1 if read is correct
*/
template<>
inline int ReadImage( const char * path, Image<RGBAColor> * im )
{
  std::vector<unsigned char> ptr;
  int w, h, depth;
  const int res = ReadImage( path, &ptr, &w, &h, &depth );
  if ( depth != 4 )
  {
    return 0;
  }
  if ( res == 1 )
  {
    RGBAColor * ptrCol = reinterpret_cast<RGBAColor*>( &ptr[0] );
    //convert raw array to Image
    ( *im ) = Eigen::Map<Image<RGBAColor>::Base>( ptrCol, h, w );
  }
  return res;
}

//--------
//-- Image Writing
//--------

/// Write image to disk, support only unsigned char based type (gray, rgb, rgba)
/**
* @brief Generic Write image
* @param filename Output image file path
* @param im Input image
* @retval 0 if there was an error during write operation
* @return non null if there was no error during write operation
*/
template<typename T>
int WriteImage( const char * filename, const Image<T>& im )
{
  const unsigned char * ptr = ( unsigned char* )( im.GetMat().data() );
  const int depth = sizeof( T ) / sizeof( unsigned char );
  std::vector<unsigned char> array( ptr , ptr + im.Width()*im.Height()*depth );
  const int w = im.Width(), h = im.Height();
  return WriteImage( filename, array, w, h, depth );
}

/**
* @brief Generic JPEG write
* @param filename Output image file path
* @param im Input image
* @param quality Output JPEG quality
* @retval 0 if there was an error during write operation
* @return non null if there was no error during write operation
* @note quality is between 0 to 100
*/
template<typename T>
int WriteJpg( const char * filename, const Image<T>& im, int quality )
{
  const unsigned char * ptr = ( unsigned char* )( im.GetMat().data() );
  const int w = im.Width(), h = im.Height();
  const int depth = sizeof( T ) / sizeof( unsigned char );
  std::vector<unsigned char> array( ptr , ptr + w * h * depth );
  return WriteJpg( filename, array, w, h, depth, quality );
}

}  // namespace image
}  // namespace openMVG

#endif  // OPENMVG_IMAGE_IMAGE_IMAGE_IO_HPP
