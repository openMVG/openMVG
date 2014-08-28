*******************
image
*******************

This module provides generic algorithms for image related tasks:

Container 
=============

An openMVG ``Image<T>`` class is a generic image container:

* an Eigen matrix aligned in row major with a template pixel type T.

  * ``template <typename T> class Image : public Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>``

Single channel image:

* a 8-bit gray image: ``Image<unsigned char>``
* a 32-bit gray image: ``Image<double>``

Multichannel image: (In order to ease usage for color images some types are already available to handle RGB and RGBA images)

* a 8-bit RGB image: ``Image<RGBColor> <=> Image<Rgb<unsigned char> >``
* a 8-bit RGBA image: ``Image<RGBAColor> <=> Image<Rgba<unsigned char> >``
* a 32-bit RGB image: ``Image<Rgb<double> >``

Image I/O 
=============

Loading and writing of 8 bits (gray and color) images are supported in the following formats:

* ppm/pgm,
* jpeg,
* png,
* tiff.

Drawing operations
===================

The following operations are available:

* lines,
* circles,
* ellipses.

An examples from openMVG/images/image_drawing_test.cpp:

.. code-block:: c++ 

  Image<unsigned char> image(10,10);
  image.fill(0);
  
  // Pixel access is done as matrix (row, line)
  int row = 2;
  int column = 4;
  image(row, column) = 127;

  // Horizontal scanline
  DrawLine( 0, 5, w-1, 5, 255, &image);

  // Circle of radius 3 and center (5,5)
  const int radius = 3;
  const int x = 5, y = 5;
  DrawCircle(x, y, radius, (unsigned char)255, &image);
        
  // Ellipse of center (5,5) and (3,0
  const int radius1 = 3, radius2 = 1, angle = 0;
  const int x = 5, y = 5;

  DrawEllipse(x, y, radius1, radius2, (unsigned char)255, &image, (double)angle);
  
  // Example with a RGB image
  Image<RGBColor> imageRGB(10,10);
  DrawCircle(x, y, radius, RGBColor(255,0,0), &imageRGB);

  

