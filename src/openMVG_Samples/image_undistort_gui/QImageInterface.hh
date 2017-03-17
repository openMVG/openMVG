#ifndef _OPENMVG_SAMPLE_IMAGE_UNDISTORT_GUI_QIMAGE_INTERFACE_HH_
#define _OPENMVG_SAMPLE_IMAGE_UNDISTORT_GUI_QIMAGE_INTERFACE_HH_

#include "openMVG/image/image.hpp"

#include <QImage>

namespace image_undistort_gui
{

/**
  * @brief Build a QImage from an openMVG (RGB) image 
  * @param img Input image 
  * @return a QImage corresponding to the input image 
  * @note this makes a deep copy of the image 
  */
QImage openMVGImageToQImage( const openMVG::image::Image<openMVG::image::RGBColor> &img );

/**
  * @brief Build a QImage from an openMVG (grayscale) image 
  * @param img Input image 
  * @return a QImage corresponding to the input image 
  * @note this makes a deep copy of the image 
  */
QImage openMVGImageToQImage( const openMVG::image::Image<unsigned char> &img );

/**
  * @brief Convert a QImage to an openMVG image 
  * @param img Input image 
  * @return openMVG image corresponding to this image 
  * @note this makes a deep copy 
  */
openMVG::image::Image<openMVG::image::RGBColor> QImageToOpenMVGImage( const QImage &img );

} // namespace image_undistort_gui

#endif