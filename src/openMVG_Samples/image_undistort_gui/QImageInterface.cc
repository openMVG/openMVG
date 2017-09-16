// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Romuald Perrot.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "QImageInterface.hh"

#include <algorithm>

#include <QRgb>
#include <QVector>

namespace image_undistort_gui
{

/**
  * @brief Build a QImage from an openMVG (RGB) image
  * @param img Input image
  * @return a QImage corresponding to the input image
  * @note this makes a deep copy of the image
  */
QImage openMVGImageToQImage( const openMVG::image::Image<openMVG::image::RGBColor> &img )
{
  const int bytePerRow = img.Width() * sizeof( openMVG::image::RGBColor );
  QImage tmp( reinterpret_cast<const unsigned char *>( img.data() ), img.Width(), img.Height(), bytePerRow, QImage::Format_RGB888 );

  // Make a deep copy
  return tmp.copy();
}

/**
  * @brief Build a QImage from an openMVG (grayscale) image
  * @param img Input image
  * @return a QImage corresponding to the input image
  * @note this makes a deep copy of the image
  */
QImage openMVGImageToQImage( const openMVG::image::Image<unsigned char> &img )
{
  QImage tmp( reinterpret_cast<const unsigned char *>( img.data() ), img.Width(), img.Height(), QImage::Format_Indexed8 );
  QVector<QRgb> colors( 256 );
  for ( int i = 0; i < 256; ++i )
  {
    colors[ i ] = qRgb(i, i, i);
  }
  tmp.setColorCount( 256 );
  tmp.setColorTable( colors );

  // Make a deep copy
  return std::move(tmp);
}

/**
  * @brief Convert a QImage to an openMVG image
  * @param img Input image
  * @return openMVG image corresponding to this image
  * @note this makes a deep copy
  */
openMVG::image::Image<openMVG::image::RGBColor> QImageToOpenMVGImage( const QImage &img )
{
  openMVG::image::Image<openMVG::image::RGBColor> res( img.width(), img.height() );

  const QImage tmp = img.convertToFormat( QImage::Format_RGB888 );

  // Make a deep copy
  std::copy( tmp.bits(), tmp.bits() + tmp.byteCount(), reinterpret_cast<unsigned char *>( res.data() ) );

  return res;
}

} // namespace image_undistort_gui
