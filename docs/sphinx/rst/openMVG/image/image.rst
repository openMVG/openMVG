*******************
image
*******************

Image Container 
===============

OpenMVG :class:`Image\<T\>` class is a generic image container based on an Eigen aligned row-major matrix of template pixel type `T`. Images can store grayscale, RGB, RGBA or custom data.

:class:`Image\<T\>` provide basic pixel read and write operation.

See examples from openMVG/images/image_test.cpp:

  .. code-block:: c++

    // A 8-bit gray image: 
    Image<unsigned char> grayscale_image_8bit;

    // A 32-bit gray image: 
    Image<double> grayscale_image_32bit;

    // Multichannel image: (use pre-defined pixel type)

    // A 8-bit RGB image: 
    Image<RGBColor> rgb_image_8bit;
    Image<Rgb<unsigned char> > rgb_image2_8bit;

    // 8-bit RGBA image
    Image<RGBAColor> rgba_image_8bit;
    Image<Rgba<unsigned char> > rgba_image2_8bit;

    // 32 bit RGB image:
    Image<Rgb<double> > rgb_image_32bit;

Image I/O 
=============

Loading and writing of 8 bits (gray and color) images are supported in the following formats:

* ppm/pgm,
* jpeg,
* png,
* tiff.

See examples from openMVG/images/image_IO_test.cpp:

  .. code-block:: c++

    // Read a grayscale image (if conversion need, it is done on the fly)
    Image<unsigned char> gray_image;
    bool bRet = ReadImage("Foo.imgExtension", &gray_image);

    // Read a color image
    Image<RGBColor> rgb_image_gray;
    bool bRet = ReadImage("Foo.imgExtension", &rgb_image);

Drawing operations
===================

The following operations are available:

* lines,
* circles,
* ellipses.

See examples from openMVG/images/image_drawing_test.cpp:

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

  

